"""Isolated execution of the published classic-SWAT fork on user inputs."""
from pathlib import Path
import calendar, hashlib, math, queue, re, shutil, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from . import policy
from .metrics import Trajectory

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def number(token):
    if isinstance(token,bytes):token=token.decode('ascii')
    token=str(token).strip()
    if not token or '*' in token:raise ValueError('Invalid Fortran number: '+repr(token))
    try:value=float(token.replace('D','E').replace('d','e'))
    except ValueError:
        match=re.fullmatch(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)',token)
        if not match:raise ValueError('Invalid Fortran number: '+repr(token)) from None
        value=float(match[1]+'E'+match[2])
    if not math.isfinite(value):raise ValueError('Nonfinite Fortran number')
    return value

def geometry(tio):
    areas={};hrus={};pattern=r'Watershed HRU:\s*(\d+)\s+Subbasin:\s*(\d+)\s+HRU:\s*(\d+)'
    for p in sorted(tio.glob('*.sub')):
        if not re.fullmatch(r'\d{9}',p.stem):continue
        lines=p.read_text(encoding='ascii',errors='replace').splitlines()
        match=re.search(r'Subbasin:\s*(\d+)',lines[0])
        if match is None:raise ValueError('Unrecognized subbasin header: '+p.name)
        area=number(next(x.split('|')[0] for x in lines if re.search(r'\|\s*SUB_KM\s*:',x)))
        sub=int(match[1])
        if sub in areas or area<=0:raise ValueError('Invalid/duplicate subbasin area')
        areas[sub]=area
    for p in sorted(tio.glob('*.hru')):
        if not re.fullmatch(r'\d{9}',p.stem):continue
        lines=p.read_text(encoding='ascii',errors='replace').splitlines();m=re.search(pattern,lines[0])
        if m is None:raise ValueError('Unrecognized HRU header: '+p.name)
        global_hru,sub,local=map(int,m.groups())
        frac=number(next(x.split('|')[0] for x in lines if re.search(r'\|\s*HRU_FR\s*:',x)))
        management=p.with_suffix('.mgt').read_text(encoding='ascii',errors='replace').splitlines()
        mm=re.search(pattern,management[0])
        if mm is None or tuple(map(int,mm.groups()))!=(global_hru,sub,local):raise ValueError('HRU/management mismatch')
        start=next(i for i,x in enumerate(management) if x.strip()=='Operation Schedule:')+1
        ids=[]
        for line in management[start:]:
            if not line.strip():continue
            op=int(line[16:18].strip() or '0');fert=int(line[19:23].strip() or '0')
            if op==3:ids.append(fert)
        agricultural=any(f in {4,44,45,47,51} for f in ids)
        if agricultural:
            if ids.count(4)!=3:raise ValueError('Each agricultural HRU needs exactly three urea operations (ID4)')
            if frac<=0 or sub not in areas or (sub,local) in hrus:raise ValueError('Invalid HRU geometry')
            hrus[sub,local]=areas[sub]*frac*100 # hectares
    if not hrus:raise ValueError('No controlled agricultural HRUs')
    return hrus

def read_actions(path,years,subs):
    with path.open(encoding='ascii') as f:
        header=f.readline().split();rows=[]
        for line in f:
            if not line.strip():continue
            values=line.split()
            if len(values)!=len(header):raise ValueError('Action column count mismatch')
            row=dict(zip(header,map(number,values)))
            if int(row['year']) in years and int(row['subbasin']) in subs:rows.append(row)
    index={(int(r['year']),int(r['subbasin'])):r for r in rows}
    if len(index)!=len(rows) or set(index)!={(int(y),int(s)) for y in years for s in subs}:
        raise ValueError('Incomplete or duplicate annual actions')
    base=[];areas=[];phase=[];oob=0
    for s in subs:
        row=index[int(years[0]),int(s)]
        base.append([row[k] for k in ['baseline_basal_N','baseline_topdress1_N','baseline_topdress2_N']]);areas.append(row['ag_area_km2']*100)
    for y in years:
        out=[]
        for j,s in enumerate(subs):
            row=index[int(y),int(s)]
            if int(row['state_year'])!=int(y) or row['manure_scale']!=1:raise ValueError('Wrong feedback timing or manure changed')
            np.testing.assert_allclose([row[k] for k in ['baseline_basal_N','baseline_topdress1_N','baseline_topdress2_N']],base[j],rtol=1e-10)
            if not math.isclose(row['ag_area_km2']*100,areas[j],rel_tol=1e-10):raise ValueError('Area changed')
            out.append([row[k] for k in ['basal_N','topdress1_N','topdress2_N']]);oob+=int(row['state_oob_count'])
        phase.append(out)
    return np.array(base),np.array(areas),np.array(phase),oob

def read_months(path,years,outlet):
    loads=[];volumes=[]
    with path.open(encoding='ascii',errors='strict') as f:
        for line in f:
            if not line.startswith('REACH'):continue
            v=[number(x) for x in line.split()[1:]]
            if len(v)!=50:raise ValueError('Expected native Rev.664 monthly reach columns')
            if int(v[0])!=outlet or not v[2].is_integer() or not 1<=v[2]<=12:continue
            i=len(loads);month=i%12+1
            if int(v[2])!=month or i>=len(years)*12:raise ValueError('Out-of-sequence monthly outlet rows')
            flow,mass=v[5],v[46]
            if flow<0 or mass<0:raise ValueError('Negative outlet flow or TN')
            loads.append(mass);volumes.append(flow*calendar.monthrange(int(years[i//12]),month)[1]*86400)
    if len(loads)!=len(years)*12:raise ValueError('Incomplete monthly reach coverage; verify IPRINT and NYSKIP')
    return np.array(loads).reshape(-1,12),np.array(volumes).reshape(-1,12)

def read_grain(path,years,subs,hrus):
    result=np.zeros((len(years),len(subs)));seen=set();expected={(int(y),s,h) for y in years for s,h in hrus}
    with path.open(encoding='ascii',errors='strict') as f:
        for line in f:
            if 'HARV/KILL' not in line:continue
            values=line.split()
            if len(values)<16 or values[7]!='HARV/KILL' or not values[0].isdigit():continue
            s,h,y=map(int,values[:3]);key=y,s,h
            if y not in years or (s,h) not in hrus:continue
            if key in seen:raise ValueError('Duplicate agricultural harvest; this release expects one annual maize harvest per HRU')
            seen.add(key);v=number(values[15])
            if v<0:raise ValueError('Negative grain harvest')
            result[int(y-years[0]),int(np.where(subs==s)[0][0])]+=v*hrus[s,h]/1000
    if seen!=expected:raise ValueError('Incomplete annual agricultural HARV/KILL records')
    return result

class SwatCase:
    name='swat'
    def __init__(self,config,workspace):
        self.config=config;self.root=Path(workspace);self.root.mkdir(parents=True,exist_ok=True)
        self.template=Path(config['swat']['template_dir']).resolve();self.exe=Path(config['swat']['executable']).resolve()
        if not self.template.is_dir() or not self.exe.is_file():raise ValueError('Provide a complete user TxtInOut and compiled executable')
        self.args=list(map(str,config['swat']['calibration_parameters']))
        if len(self.args)!=27 or not all(math.isfinite(float(x)) for x in self.args):raise ValueError('Exactly27 finite calibration values required')
        self.output_id=str(config['swat'].get('output_id','959999'))
        if not re.fullmatch(r'\d{1,6}',self.output_id):raise ValueError('Numeric output ID up to6 digits required')
        self.years=np.arange(config['evaluation_years'][0],config['evaluation_years'][1]+1)
        self.hrus=geometry(self.template);self.subs=np.array(sorted({s for s,h in self.hrus}))
        self.inputs={p.name:sha(p) for p in self.template.iterdir() if p.is_file() and self._is_input(p)}
        if 'file.cio' not in self.inputs:raise ValueError('file.cio missing')
        self._check_cio();self.slots=queue.Queue();self.slot_paths=[]
        self.engine_hash=sha(self.exe)
        for i in range(int(config.get('workers',1))):
            tio=self.root/('worker_%02d'%i)/'TxtInOut'
            if tio.exists():raise ValueError('Worker workspace exists; optimizer uses a fresh worker root on resume')
            tio.mkdir(parents=True)
            for name in self.inputs:shutil.copy2(self.template/name,tio/name)
            # The inherited Windows source routes suppressed HRU output to NUL.
            # On POSIX this must be a null-device link, not a large disk file.
            if not hasattr(subprocess,'CREATE_NO_WINDOW'):(tio/'NUL').symlink_to('/dev/null')
            self.slots.put(tio);self.slot_paths.append(tio)

    @staticmethod
    def _is_input(path):
        n=path.name.lower()
        generated=re.fullmatch(r'\d{1,6}(?:input\.std|watout\.dat|(?:hru|rch|rsv|sub)\.dat|parameter_audit\.txt|chan\.deg|swat\.qst)',n)
        return not (generated or 'output' in n or n.startswith(('dps_','calibration_')) or n in {'fin.fin','run.log','nul'} or
                    path.suffix.lower() in {'.exe','.dll','.obj','.mod','.o','.log','.out','.pkl','.json','.tsv'})

    def _check_cio(self):
        fields={}
        for line in (self.template/'file.cio').read_text(encoding='ascii',errors='replace').splitlines():
            if '|' in line:
                value,label=line.split('|',1);m=re.match(r'\s*([A-Z0-9_]+)\s*:',label)
                if m:fields[m[1]]=value.strip()
        start,end=self.config['evaluation_years'];warm=self.config['warmup_years']
        expected={'IYR':warm[0],'NBYR':end-warm[0]+1,'NYSKIP':start-warm[0],'IPRINT':0,'IDAF':1,
                  'IDAL':366 if calendar.isleap(end) else 365}
        for key,value in expected.items():
            if int(float(fields.get(key,'nan')))!=value:raise ValueError('file.cio '+key+' must equal '+str(value))
        if warm[1]!=start-1:raise ValueError('Warmup/evaluation mismatch')

    def fingerprint(self):
        return {'backend':'swat','inputs':self.inputs,'executable_sha256':self.engine_hash,'configuration':self.config}

    def simulate(self,theta,baseline=False):
        tio=self.slots.get()
        try:
            for name,digest in self.inputs.items():
                if sha(tio/name)!=digest:raise ValueError('Worker input changed: '+name)
            for name in ['fin.fin',self.output_id+'output.rch',self.output_id+'output.mgt','dps_actions.tsv']:
                p=tio/name
                if p.is_file():p.unlink()
            c=self.config
            policy.write(tio/'dps_policy.in',theta,c['state_lower'],c['state_upper'],mode=0 if baseline else 1,
                         start=c['evaluation_years'][0],annual_bounds=c['annual_bounds'],phase_lower=c['phase_lower'],
                         phase_upper=c['phase_upper'],phase_reference=c['phase_reference'],outlet=int(c['outlet_reach']))
            kwargs={'creationflags':subprocess.CREATE_NO_WINDOW} if hasattr(subprocess,'CREATE_NO_WINDOW') else {}
            with (tio.parent/'run.log').open('w',encoding='utf-8') as log:
                result=subprocess.run([str(self.exe),*self.args,self.output_id],cwd=tio,stdout=log,stderr=subprocess.STDOUT,
                    timeout=float(c['swat'].get('timeout_seconds',900)),**kwargs)
            if result.returncode or not (tio/'fin.fin').exists():raise RuntimeError('SWAT failed; see '+str(tio.parent/'run.log'))
            base,areas,phase,oob=read_actions(tio/'dps_actions.tsv',self.years,self.subs)
            for j,s in enumerate(self.subs):
                expected=sum(v for (sub,h),v in self.hrus.items() if sub==s)
                if not math.isclose(areas[j],expected,rel_tol=2e-6):raise ValueError('Input/output area mismatch')
            mass,volume=read_months(tio/(self.output_id+'output.rch'),self.years,int(c['outlet_reach']))
            grain=read_grain(tio/(self.output_id+'output.mgt'),self.years,self.subs,self.hrus)
            for name,digest in self.inputs.items():
                if sha(tio/name)!=digest:raise ValueError('SWAT input changed during execution: '+name)
            return Trajectory(self.years,self.subs,areas,base,phase,grain,mass,volume,oob).validate()
        finally:self.slots.put(tio)
