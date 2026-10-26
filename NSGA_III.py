"""Eight-state fertilizer DPS + NSGA-III using a real classic-SWAT executable.

1. Install requirements.txt and a Fortran compiler.
2. python NSGA_III.py --build-only
3. python NSGA_III.py --population 8 --generations 2 --workers 2

CASE below configures the supplied two-subbasin TxtInOut example.
No surrogate model or research observations are used.
"""
from pathlib import Path
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,calendar,csv,hashlib,json,math,os,platform,queue,random,re,shutil,signal,subprocess,sys,threading,time
import dill as pickle
import numpy as np
import pymoo
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.core.problem import Problem
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from pymoo.indicators.hv import HV
from pymoo.termination import get_termination

ROOT=Path(__file__).resolve().parent
POPULATION=100
GENERATIONS=1400
WORKERS=8
SEED=20261007

# Illustrative fixed inputs, not fitted study-basin values. To use your own
# project, update years, outlet, fixed state bounds and hydrological values.
CASE={'backend': 'swat',
 'evaluation_years': [2002, 2005],
 'warmup_years': [2000, 2001],
 'outlet_reach': 2,
 'state_lower': [0, 0, 0, 0, 0, 0, 0, 0],
 'state_upper': [500, 10, 300, 600, 20000, 80, 1.5, 600],
 'annual_bounds': [0.2, 1.2],
 'phase_lower': [0.3, 0.1, 0.3],
 'phase_upper': [0.5, 0.3, 0.5],
 'phase_reference': [0.4, 0.2, 0.4],
 'grain_floor_fraction': 0.95,
 'TN_limit_mg_L': 4.0,
 'TN_max_exceedance_fraction': 0.25,
 'HV_reference_scaled': [1.2, 0, 1.2, 0],
 'swat': {'calibration_parameters': [2,
                                     0,
                                     0.04,
                                     0.5,
                                     0,
                                     0,
                                     0,
                                     4,
                                     1,
                                     0.05,
                                     0.95,
                                     31,
                                     1,
                                     -6.5,
                                     0.5,
                                     1,
                                     4,
                                     0.02,
                                     100,
                                     0.005,
                                     4,
                                     0.5,
                                     100,
                                     0.15,
                                     4,
                                     1,
                                     1],
          'output_id': '1',
          'timeout_seconds': 120}}


# ---- RBF policy and fertilizer allocation ----
B,K,N,P=8,3,12,231

SCHEMA='DPS_SUBBASIN_SPLITN_V5'

STATES=('soil_nitrate_N','soil_ammonium_N','soil_fresh_org_N',
        'soil_active_org_N','soil_stable_org_N','aquifer_NO3_N',
        'soil_water','shallow_aquifer_water')

PHASES=('basal','topdress1','topdress2')

REGIONAL_SHARES=np.array([.4,.2,.4])

PHASE_LO=np.array([.3,.1,.3])

PHASE_HI=np.array([.5,.3,.5])

def policy_vector(bias=None,weights=None,centers=None,widths=None):
    if bias is None:bias=REGIONAL_SHARES/1.2
    return np.r_[np.asarray(bias).ravel(),
        np.zeros((N,K)).ravel() if weights is None else np.asarray(weights).ravel(),
        np.full((N,B),.5).ravel() if centers is None else np.asarray(centers).ravel(),
        np.ones((N,B)).ravel() if widths is None else np.asarray(widths).ravel()]

def unpack_policy(theta):
    t=np.asarray(theta,dtype=float)
    if t.shape!=(P,) or not np.isfinite(t).all():raise ValueError('Expected 231 finite policy parameters')
    a=t[:K];w=t[K:K+N*K].reshape(N,K)
    c=t[K+N*K:K+N*K+N*B].reshape(N,B);r=t[K+N*K+N*B:].reshape(N,B)
    if not ((a>=0)&(a<=1)).all() or not (abs(w)<=1).all():raise ValueError('Bias/weight range')
    if not ((c>=0)&(c<=1)).all() or not ((r>=.01)&(r<=1)).all():raise ValueError('Center/width range')
    return a,w,c,r

def validate_phase_bounds(lower,upper,reference):
    lower=np.asarray(lower,dtype=float);upper=np.asarray(upper,dtype=float);reference=np.asarray(reference,dtype=float)
    if lower.shape!=upper.shape or lower.shape!=reference.shape or lower.shape!=(3,):raise ValueError('Three phase bounds required')
    if not np.isfinite(np.r_[lower,upper,reference]).all():raise ValueError('Nonfinite phase bounds')
    if not ((lower>=0)&(upper<=1)&(lower<=upper)).all() or lower.sum()>1 or upper.sum()<1:raise ValueError('Infeasible phase bounds')
    if abs(reference.sum()-1)>1e-12 or not ((reference>=lower)&(reference<=upper)).all():raise ValueError('Reference outside phase space')
    return lower,upper,reference

def project_shares(shares,lower=PHASE_LO,upper=PHASE_HI):
    """Independent active-set enumeration; not Fortran's bisection algorithm."""
    import itertools
    shares=np.asarray(shares,dtype=float);lower=np.asarray(lower);upper=np.asarray(upper)
    if shares.shape!=(3,) or not np.isfinite(shares).all():raise ValueError('Invalid shares')
    if ((shares>=lower)&(shares<=upper)).all() and abs(shares.sum()-1)<1e-12:return shares.copy()
    best=None;score=float('inf')
    # Status -1: lower bound, 0: free, 1: upper bound; 27 possible active sets.
    for status in itertools.product((-1,0,1),repeat=3):
        free=[i for i in range(3) if status[i]==0]
        fixed=[i for i in range(3) if status[i]!=0]
        out=np.zeros(3)
        for i in fixed:out[i]=lower[i] if status[i]==-1 else upper[i]
        if free:
            shift=(sum(shares[i] for i in free)+sum(out[i] for i in fixed)-1)/len(free)
            for i in free:out[i]=shares[i]-shift
        elif abs(out.sum()-1)>1e-12:continue
        if (out<lower-1e-12).any() or (out>upper+1e-12).any() or abs(out.sum()-1)>1e-12:continue
        value=float(np.sum((out-shares)**2))
        if value<score:best=out;score=value
    if best is None:raise ValueError('Share projection infeasible')
    return best

def allocate(u,baseline_phase_N,annual_bounds=(.2,1.2),phase_lower=PHASE_LO,phase_upper=PHASE_HI,phase_reference=REGIONAL_SHARES):
    u=np.asarray(u,dtype=float);base=np.asarray(baseline_phase_N,dtype=float)
    lower,upper=annual_bounds
    if u.shape!=(3,) or base.shape!=(3,) or not np.isfinite(np.r_[u,base,lower,upper]).all():raise ValueError('Budget inputs')
    if not ((u>=0)&(u<=1)).all() or not (base>0).all():raise ValueError('Budget ranges')
    if not 0<=lower<=1<=upper or not lower<upper:raise ValueError('Annual bounds must include baseline')
    # The reference is the sum of existing chemical operation doses, never manure.
    ref=float(base.sum());sum_u=float(u.sum());trial=upper*sum_u
    scale=float(np.clip(trial,lower,upper));total=ref*scale
    fallback=sum_u==0
    phase_lower,phase_upper,phase_reference=validate_phase_bounds(phase_lower,phase_upper,phase_reference)
    trial_shares=phase_reference.copy() if fallback else u/sum_u
    shares=project_shares(trial_shares,phase_lower,phase_upper)
    target=total*shares
    return dict(target=target,phase_scale=target/base,annual_scale=scale,
                total=float(target.sum()),trial_annual_scale=trial,
                budget_clip=-1 if trial<lower else (1 if trial>upper else 0),zero_fallback=int(fallback),
                trial_shares=trial_shares,shares=shares,
                share_repaired=int((trial_shares<phase_lower).any() or (trial_shares>phase_upper).any()))

def evaluate_policy(theta,state,lo,hi,baseline_phase_N,annual_bounds=(.2,1.2),phase_lower=PHASE_LO,phase_upper=PHASE_HI,phase_reference=REGIONAL_SHARES):
    a,w,c,r=unpack_policy(theta)
    state=np.asarray(state);lo=np.asarray(lo);hi=np.asarray(hi)
    if state.shape!=lo.shape or lo.shape!=hi.shape or state.shape!=(B,):raise ValueError('Eight state scales required')
    if not np.isfinite(np.r_[state,lo,hi]).all() or not (hi>lo).all():raise ValueError('Invalid state scales')
    z=np.clip((state-lo)/(hi-lo),0,1)
    phi=np.exp(-np.sum(((z-c)/r)**2,axis=1))
    raw=a+phi@w;u=np.clip(raw,0,1)
    result=allocate(u,baseline_phase_N,annual_bounds,phase_lower,phase_upper,phase_reference)
    result.update(z=z,phi=phi,raw=raw,u=u)
    return result

def write_policy(path,theta,lo,hi,mode=1,start=2008,annual_bounds=(.2,1.2),phase_lower=PHASE_LO,phase_upper=PHASE_HI,phase_reference=REGIONAL_SHARES,outlet=21):
    a,w,c,r=unpack_policy(theta)
    evaluate_policy(theta,lo,lo,hi,np.ones(3),annual_bounds,phase_lower,phase_upper,phase_reference)
    if mode not in (0,1):raise ValueError('Only baseline mode 0 and RBF mode 1 supported')
    if not isinstance(outlet,int) or outlet<1:raise ValueError('Positive integer outlet required')
    rows=[SCHEMA,f'{mode} {start} {outlet}',f'{B} {K} {N} {P}',
        ' '.join(f'{x:.17g}' for x in annual_bounds),
        ' '.join(f'{x:.17g}' for x in phase_lower),
        ' '.join(f'{x:.17g}' for x in phase_upper),
        ' '.join(f'{x:.17g}' for x in phase_reference),
        ' '.join(f'{x:.17g}' for x in lo),' '.join(f'{x:.17g}' for x in hi),
        ' '.join(f'{x:.17g}' for x in a)]
    rows += [' '.join(f'{x:.17g}' for x in row) for block in (w,c,r) for row in block]
    Path(path).write_text('\n'.join(rows)+'\n',encoding='ascii')

def read_policy(path):
    rows=Path(path).read_text(encoding='ascii').splitlines()
    if len(rows)!=46 or rows[0]!=SCHEMA:raise ValueError('Expected SPLITN V5, 46 rows')
    if tuple(map(int,rows[2].split()))!=(B,K,N,P):raise ValueError('Policy dimensions')
    mode,start,outlet=map(int,rows[1].split())
    expected=[2,3,3,3,8,8,3]+[3]*N+[8]*(2*N)
    if any(len(row.split())!=width for row,width in zip(rows[3:],expected)):raise ValueError('Policy row width')
    bounds,phase_lower,phase_upper,phase_reference,lo,hi=[np.array(list(map(float,row.split()))) for row in rows[3:9]]
    theta=np.array([float(x) for row in rows[9:] for x in row.split()]);unpack_policy(theta)
    evaluate_policy(theta,lo,lo,hi,np.ones(3),bounds,phase_lower,phase_upper,phase_reference)
    if mode not in (0,1) or outlet<1:raise ValueError('Policy mode/outlet')
    return dict(mode=mode,start=start,outlet=outlet,bounds=bounds,lo=lo,hi=hi,theta=theta,
                phase_lower=phase_lower,phase_upper=phase_upper,phase_reference=phase_reference)

# ---- Objectives and constraints ----
@dataclass
class Trajectory:
    years: np.ndarray
    subbasins: np.ndarray
    area_ha: np.ndarray
    baseline_phase_N: np.ndarray  # [subbasin, phase], kg N/ha
    phase_N: np.ndarray           # [year, subbasin, phase], kg N/ha
    grain_t: np.ndarray           # [year, subbasin], dry tonnes
    monthly_TN_kg: np.ndarray     # [year, month], kg N
    monthly_volume_m3: np.ndarray # [year, month], m3
    state_oob_count: int = 0

    def validate(self):
        y=np.asarray(self.years);s=np.asarray(self.subbasins);ny,ns=len(y),len(s)
        if ny<2 or ns<1 or not np.array_equal(y,np.arange(y[0],y[-1]+1)):
            raise ValueError('At least two consecutive evaluation years required')
        if len(set(s))!=ns:raise ValueError('Duplicate subbasin')
        expected={'area_ha':(ns,), 'baseline_phase_N':(ns,3), 'phase_N':(ny,ns,3),
                  'grain_t':(ny,ns), 'monthly_TN_kg':(ny,12), 'monthly_volume_m3':(ny,12)}
        for name,shape in expected.items():
            arr=np.asarray(getattr(self,name),dtype=float)
            if arr.shape!=shape or not np.isfinite(arr).all() or (arr<0).any():
                raise ValueError('Invalid trajectory field: '+name)
            setattr(self,name,arr)
        if (self.area_ha<=0).any() or (self.baseline_phase_N<=0).any():
            raise ValueError('Positive agricultural area and baseline event doses required')
        if np.any((self.monthly_volume_m3==0)&(self.monthly_TN_kg>0)):
            raise ValueError('Positive TN mass with zero outlet water volume')
        return self

def evaluate_metrics(current,baseline,config):
    """J=[TN,-grain,executed-action smoothness,-worst local retention]."""
    current.validate();baseline.validate()
    for name in ['years','subbasins','area_ha','baseline_phase_N']:
        if not np.array_equal(getattr(current,name),getattr(baseline,name)):
            raise ValueError('Policy/baseline comparison has mismatched '+name)
    base_local=baseline.grain_t.sum(axis=0)
    if (base_local<=0).any():raise ValueError('Positive baseline grain in every managed subbasin required')
    retention=current.grain_t.sum(axis=0)/base_local
    grain=float(current.grain_t.sum()/len(current.years))
    base_grain=float(baseline.grain_t.sum()/len(baseline.years))
    ratio=grain/base_grain
    # Normalize each executed phase change by the fixed LOCAL ANNUAL baseline,
    # sum the 3 phases, area-weight subbasins, average T-1 transitions, x100.
    change=np.abs(np.diff(current.phase_N,axis=0))/current.baseline_phase_N.sum(axis=1)[None,:,None]
    smooth=float(np.mean(np.sum(change.sum(axis=2)*current.area_ha[None,:],axis=1)/current.area_ha.sum())*100)
    concentration=np.divide(current.monthly_TN_kg*1000,current.monthly_volume_m3,
                            out=np.zeros_like(current.monthly_TN_kg),where=current.monthly_volume_m3>0)
    limit=float(config.get('TN_limit_mg_L',4.0));floor=float(config.get('grain_floor_fraction',.95))
    allowed=float(config.get('TN_max_exceedance_fraction',.25))
    if not (limit>0 and 0<floor<=1 and 0<allowed<=1):raise ValueError('Invalid constraint settings')
    count=int(np.count_nonzero(concentration>limit));months=int(concentration.size)
    # Strict >95% and strict <25%. For 516 months, max_count=128, not129.
    maximum=math.ceil(allowed*months)-1
    g=[(float(np.nextafter(floor,math.inf))-ratio)/floor,(count-maximum)/months]
    load=float(current.monthly_TN_kg.sum()/len(current.years)/1000)
    base_load=float(baseline.monthly_TN_kg.sum()/len(baseline.years)/1000)
    if base_load<=0:raise ValueError('Positive paired baseline TN required')
    f=[load,-grain,smooth,-float(retention.min())]
    scaled=[load/base_load,-ratio,smooth/100,-float(retention.min())]
    chem=float(np.mean(np.sum(current.phase_N.sum(axis=2)*current.area_ha[None,:],axis=1))/1000)
    return {'F':f,'F_scaled':scaled,'G':g,'feasible':bool(all(x<=0 for x in g)),
        'TN_t_per_year':load,'grain_t_dry_per_year':grain,'smoothness_pp_per_year':smooth,
        'grain_retention':ratio,'worst_subbasin_retention':float(retention.min()),
        'worst_subbasin':int(current.subbasins[int(np.argmin(retention))]),
        'TN_exceedance_count':count,'months':months,'TN_max_exceedance_count':maximum,
        'TN_exceedance_fraction':count/months,'chemical_N_t_per_year':chem,
        'state_oob_count':int(current.state_oob_count),'transitions':len(current.years)-1,
        'CV':float(sum(max(0,x) for x in g))}

def write_trajectory(directory,trajectory):
    """Portable output tables; no private data are bundled with the release."""
    from pathlib import Path
    import csv
    p=Path(directory);p.mkdir(parents=True,exist_ok=True);t=trajectory.validate()
    with (p/'actions.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['year','subbasin','area_ha','baseline_annual_N_kg_ha','basal_N_kg_ha','topdress1_N_kg_ha','topdress2_N_kg_ha','annual_scale'])
        for iy,y in enumerate(t.years):
            for isub,s in enumerate(t.subbasins):
                ref=t.baseline_phase_N[isub].sum();dose=t.phase_N[iy,isub]
                w.writerow([int(y),int(s),t.area_ha[isub],ref,*dose,dose.sum()/ref])
    with (p/'annual_grain.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['year','subbasin','grain_t_dry'])
        for iy,y in enumerate(t.years):
            for j,s in enumerate(t.subbasins):w.writerow([int(y),int(s),t.grain_t[iy,j]])
    with (p/'monthly_outlet.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['year','month','TN_kg','volume_m3','TN_mg_L'])
        for iy,y in enumerate(t.years):
            for m in range(12):
                mass=t.monthly_TN_kg[iy,m];vol=t.monthly_volume_m3[iy,m]
                w.writerow([int(y),m+1,mass,vol,mass*1000/vol if vol else 0.])

# ---- Complete SWAT execution ----
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
            # Native mean-annual summary MON is decimal (e.g.4.0), not a month.
            if int(v[0])!=outlet or not line.split()[3].isdigit() or not 1<=v[2]<=12:continue
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
            write_policy(tio/'dps_policy.in',theta,c['state_lower'],c['state_upper'],mode=0 if baseline else 1,
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

# ---- NSGA-III search and checkpointing ----
def atomic_json(path,value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp.%s.%s'%(os.getpid(),threading.get_ident()))
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8');os.replace(tmp,p)

def key_for(theta):return hashlib.sha256(np.asarray(theta,dtype='<f8').tobytes()).hexdigest()

def bounds():
    return (np.r_[np.zeros(K),-np.ones(N*K),np.zeros(N*B),np.full(N*B,.01)],
            np.ones(P))

def initial_population(size,seed):
    rng=np.random.default_rng(seed);lo,hi=bounds();x=rng.uniform(lo,hi,(size,P))
    # Generic constant policies, not optimized vectors from a research basin.
    for i,scale in enumerate([1.,.8,.6,.4,1.2,.2]):
        if i>=size:break
        x[i]=policy_vector(REGIONAL_SHARES*scale/1.2)
    return x

def load_config(args):
    c=json.loads(json.dumps(CASE))
    c['swat']['template_dir']=str(args.txtinout.resolve())
    c['swat']['executable']=str(args.exe.resolve())
    lo=np.asarray(c['state_lower'],dtype=float);hi=np.asarray(c['state_upper'],dtype=float)
    if lo.shape!=(B,) or hi.shape!=(B,) or not np.isfinite(np.r_[lo,hi]).all() or not (hi>lo).all():
        raise ValueError('Eight fixed finite state ranges with positive widths required')
    evaluate_policy(policy_vector(),lo,lo,hi,np.ones(3),c['annual_bounds'],c['phase_lower'],c['phase_upper'],c['phase_reference'])
    if c['annual_bounds']!=[.2,1.2] or c['phase_lower']!=[.3,.1,.3] or c['phase_upper']!=[.5,.3,.5]:
        raise ValueError('Expected manuscript annual and phase bounds')
    if (c['grain_floor_fraction'],c['TN_limit_mg_L'],c['TN_max_exceedance_fraction'])!=(.95,4.,.25):
        raise ValueError('Expected grain>95% and TN>4mg/L months<25%')
    ref=np.asarray(c['HV_reference_scaled'],dtype=float)
    if ref.shape!=(4,) or not np.isfinite(ref).all():raise ValueError('Four finite HV reference values required')
    return c

def save_trajectory(path,t):
    np.savez_compressed(path,**{name:getattr(t,name) for name in Trajectory.__dataclass_fields__})

def load_trajectory(path):
    with np.load(path,allow_pickle=False) as data:
        return Trajectory(**{name:data[name] for name in Trajectory.__dataclass_fields__}).validate()

def write_archive(output,archive):
    fields=['key','TN_t_per_year','grain_t_dry_per_year','smoothness_pp_per_year','worst_subbasin_retention',
            'grain_retention','TN_exceedance_fraction','TN_exceedance_count','months','chemical_N_t_per_year','state_oob_count']
    temp=output/'pareto.csv.tmp'
    with temp.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row in archive:w.writerow({key:row[key] for key in fields})
    os.replace(temp,output/'pareto.csv')
    np.savez_compressed(output/'pareto_policies.npz',theta=np.asarray([r['theta'] for r in archive]).reshape(-1,P),
                        keys=np.asarray([r['key'] for r in archive],dtype='U64'))

def archive_update(archive,records,X):
    merged={r['key']:r for r in archive}
    for row,theta in zip(records,X):
        if row['feasible']:merged[row['key']]=dict(row,theta=np.asarray(theta).tolist())
    rows=list(merged.values())
    if not rows:return []
    indices=NonDominatedSorting().do(np.array([r['F_scaled'] for r in rows]),only_non_dominated_front=True)
    return [rows[int(i)] for i in indices]

def run(args):
    config=load_config(args);out=Path(args.output).resolve()
    if args.population<4 or args.generations<1 or args.workers<1:raise ValueError('Invalid population/generations/workers')
    if out.exists() and any(out.iterdir()) and not args.resume:raise ValueError('Output exists; use --resume or a new --output directory')
    out.mkdir(parents=True,exist_ok=True)
    lock=out/'optimizer.lock'
    try:handle=lock.open('x',encoding='ascii')
    except FileExistsError:raise ValueError('Optimizer lock exists; inspect the PID before removing a stale lock') from None
    handle.write(str(os.getpid()));handle.close()
    try:return _run(args,config,out)
    finally:lock.unlink(missing_ok=True)

def _run(args,config,out):
    code_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}
    signature_data={'config':config,'population':args.population,'seed':args.seed,'code':code_hashes,
                    'numpy':np.__version__,'pymoo':pymoo.__version__,'dill':pickle.__version__,'python':platform.python_version()}
    signature=hashlib.sha256(json.dumps(signature_data,sort_keys=True).encode()).hexdigest()
    checkpoint=out/'checkpoint.pkl';previous=None
    if args.resume:
        if not checkpoint.exists():raise ValueError('No checkpoint to resume')
        # Checkpoints are local Python objects; only load files from a trusted run.
        with checkpoint.open('rb') as f:previous=pickle.load(f)
        if previous['signature']!=signature:raise ValueError('Configuration, source, environment, population or seed changed; start a new run')
        if previous['completed']>=args.generations:
            print('[DONE] Requested generation count already completed',flush=True);return
    config=dict(config);config['workers']=args.workers
    runtime=out/('workers_'+str(time.time_ns()))
    backend=SwatCase(config,runtime)
    model_fingerprint=backend.fingerprint()
    if previous is not None and previous['model_fingerprint']!=model_fingerprint:
        # Worker count may change without changing the scientific model.
        old=dict(previous['model_fingerprint']);new=dict(model_fingerprint)
        old['configuration']=dict(old['configuration']);new['configuration']=dict(new['configuration'])
        old['configuration'].pop('workers',None);new['configuration'].pop('workers',None)
        if old!=new:raise ValueError('Model inputs or executable changed since checkpoint')
    baseline_path=out/'baseline.npz'
    if previous is not None:
        if hashlib.sha256(baseline_path.read_bytes()).hexdigest()!=previous['baseline_sha256']:raise ValueError('Baseline changed')
        baseline=load_trajectory(baseline_path)
    else:
        print('[INIT] Running paired baseline',flush=True)
        baseline=backend.simulate(policy_vector(),baseline=True);save_trajectory(baseline_path,baseline)
        write_trajectory(out/'baseline',baseline)
    baseline_hash=hashlib.sha256(baseline_path.read_bytes()).hexdigest()
    baseline_result=evaluate_metrics(baseline,baseline,config);atomic_json(out/'baseline_metrics.json',baseline_result)
    partitions=max(h for h in range(1,args.population) if math.comb(h+3,3)<=args.population)
    refs=get_reference_directions('das-dennis',4,n_partitions=partitions)
    xl,xu=bounds();problem=Problem(n_var=P,n_obj=4,n_ieq_constr=2,xl=xl,xu=xu)
    if previous is None:
        algorithm=NSGA3(ref_dirs=refs,pop_size=args.population,sampling=initial_population(args.population,args.seed),
            crossover=SBX(prob=.9,eta=15),mutation=PM(prob=1.,prob_var=1/P,eta=20),eliminate_duplicates=True)
        algorithm.setup(problem,seed=args.seed,termination=('n_gen',args.generations),verbose=False)
        completed=0;archive=[];history=[];requests=0;model_runs=0;cache_hits=0;prior_seconds=0.
    else:
        algorithm=previous['algorithm'];algorithm.termination=get_termination('n_gen',args.generations)
        completed=previous['completed'];archive=previous['archive'];history=previous['history']
        requests=previous['requests'];model_runs=previous['model_runs'];cache_hits=previous['cache_hits'];prior_seconds=history[-1]['elapsed_seconds']
    ref=np.asarray(config['HV_reference_scaled']);indicator=HV(ref_point=ref);stop=[False]
    def request_stop(signum,frame):
        stop[0]=True;print('[STOP] Finish current generation and checkpoint',flush=True)
    old_signal=signal.signal(signal.SIGINT,request_stop)
    cache=out/'cache';cache.mkdir(exist_ok=True)
    def one(theta):
        key=key_for(theta);p=cache/(key+'.json')
        if p.exists():
            result=json.loads(p.read_text(encoding='utf-8'))
            if result.get('signature')!=signature or result.get('key')!=key:raise ValueError('Stale evaluation cache')
            return result,True
        trajectory=backend.simulate(theta);result=evaluate_metrics(trajectory,baseline,config)
        result.update(key=key,signature=signature);atomic_json(p,result)
        return result,False
    manifest={'schema':'SWAT_DPS_NSGAIII_V2','signature':signature,'backend':config['backend'],'population':args.population,
        'generations':args.generations,'seed':args.seed,'workers':args.workers,'B':B,'K':K,'N':N,'P':P,
        'objectives':4,'constraints':2,'reference_directions':len(refs),'reference_partitions':partitions,
        'crossover':{'type':'SBX','probability':.9,'eta':15},'mutation':{'type':'PM','probability':1.,'per_variable':1/P,'eta':20},
        'configuration':config,'runtime':{k:signature_data[k] for k in ['python','numpy','pymoo','dill']},
        'HV_reference_scaled':ref.tolist(),'HV_scope':'feasible cumulative nondominated archive, unchanged reference and scales',
        'generation_count_includes_initialization':True,'model':'classic SWAT2016 Rev.664'}
    atomic_json(out/'manifest.json',manifest);start=time.monotonic()
    if previous is not None:
        # pymoo0.6.1.x uses global RNGs;0.6.2 also stores an algorithm RNG.
        # Preserve both so a resumed run follows the same candidate sequence.
        np.random.set_state(previous['numpy_random_state'])
        random.setstate(previous['python_random_state'])
    print('[START] NSGA-III population=%d generations=%d workers=%d B=8 K=3 N=12 P=231'%(args.population,args.generations,args.workers),flush=True)
    try:
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            while completed<args.generations and not stop[0]:
                infills=algorithm.ask();X=np.asarray(infills.get('X'));records=[None]*len(X)
                futures={executor.submit(one,x):i for i,x in enumerate(X)}
                for f in as_completed(futures):
                    i=futures[f];records[i],hit=f.result();cache_hits+=int(hit);model_runs+=int(not hit)
                infills.set('F',np.asarray([r['F_scaled'] for r in records]));infills.set('G',np.asarray([r['G'] for r in records]))
                algorithm.tell(infills=infills);completed+=1;requests+=len(X)
                archive=archive_update(archive,records,X)
                F=np.asarray([r['F_scaled'] for r in archive]).reshape(-1,4)
                covered=F[np.all(F<ref,axis=1)];hv=float(indicator(covered)) if len(covered) else 0.
                pop_G=algorithm.pop.get('G');feasible=int(np.count_nonzero(np.all(pop_G<=0,axis=1)))
                row={'generation':completed,'evaluation_requests':requests,'model_evaluations':model_runs,'cache_hits':cache_hits,
                     'feasible_population':feasible,'archive_size':len(archive),'HV':hv,'outside_HV_reference':len(F)-len(covered),
                     'minimum_CV':float(np.maximum(pop_G,0).sum(axis=1).min()),'elapsed_seconds':prior_seconds+time.monotonic()-start}
                history.append(row)
                with (out/'history.csv').open('w',newline='',encoding='utf-8') as f:
                    w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();w.writerows(history)
                state=dict(signature=signature,model_fingerprint=model_fingerprint,baseline_sha256=baseline_hash,algorithm=algorithm,
                           completed=completed,archive=archive,history=history,requests=requests,model_runs=model_runs,cache_hits=cache_hits,
                           numpy_random_state=np.random.get_state(),python_random_state=random.getstate())
                temp=out/'checkpoint.pkl.tmp'
                with temp.open('wb') as f:pickle.dump(state,f,protocol=pickle.HIGHEST_PROTOCOL)
                if checkpoint.exists():shutil.copy2(checkpoint,out/'checkpoint.previous.pkl')
                os.replace(temp,checkpoint);write_archive(out,archive)
                atomic_json(out/'status.json',{'status':'running','completed_generations':completed,'requested_generations':args.generations,
                    'evaluation_requests':requests,'feasible_archive_size':len(archive),'HV':hv})
                print('[GEN %04d/%d] eval=%d feasible=%d/%d archive=%d HV=%.8f CV_min=%.5g elapsed=%.1fs'%
                      (completed,args.generations,requests,feasible,len(algorithm.pop),len(archive),hv,row['minimum_CV'],row['elapsed_seconds']),flush=True)
                if (out/'STOP').exists():stop[0]=True
        status='completed' if completed==args.generations else 'stopped_after_generation'
        atomic_json(out/'status.json',{'status':status,'completed_generations':completed,'requested_generations':args.generations,
                    'evaluation_requests':requests,'feasible_archive_size':len(archive),'HV':history[-1]['HV'] if history else 0.})
        if archive:
            best=min(archive,key=lambda r:r['TN_t_per_year'])
            theta=np.asarray(best['theta']);atomic_json(out/'example_policy.json',theta.tolist())
            write_policy(out/'dps_policy.in',theta,config['state_lower'],config['state_upper'],start=config['evaluation_years'][0],
                         annual_bounds=config['annual_bounds'],phase_lower=config['phase_lower'],phase_upper=config['phase_upper'],
                         phase_reference=config['phase_reference'],outlet=int(config['outlet_reach']))
            executed=backend.simulate(theta);result=evaluate_metrics(executed,baseline,config)
            np.testing.assert_allclose(result['F'],best['F'],rtol=1e-10,atol=1e-10)
            atomic_json(out/'example_policy_metrics.json',result);write_trajectory(out/'example_policy',executed)
        print('[%s] %s'%(status.upper(),out),flush=True)
    except Exception as e:
        atomic_json(out/'status.json',{'status':'failed','completed_generations':completed,'error':str(e)})
        raise
    finally:signal.signal(signal.SIGINT,old_signal)


# ---- Compiler and command-line entry ----

def build_swat(compiler=None):
    """Compile all SWAT sources; no precompiled binary is distributed."""
    compiler=compiler or ('ifx' if os.name=='nt' else 'gfortran')
    src=ROOT/'SWAT_Source';out=ROOT/'build';out.mkdir(exist_ok=True)
    exe=out/('SWAT_DPS.exe' if os.name=='nt' else 'SWAT_DPS')
    module=src/'modparm.f'
    sources=sorted(p for p in src.iterdir() if p.suffix.lower() in {'.f','.f90'} and p.name!='modparm.f')
    log_path=out/'compile.log'
    with log_path.open('w',encoding='utf-8') as log:
        if compiler=='ifx':
            if os.name!='nt':raise ValueError('Use gfortran on Linux')
            setup=Path(os.environ.get('ONEAPI_ROOT',r'C:\Program Files (x86)\Intel\oneAPI'))/'setvars.bat'
            if not setup.is_file():raise ValueError('Install Intel oneAPI Fortran or set ONEAPI_ROOT')
            (out/'sources.rsp').write_text('\n'.join('"'+str(p)+'"' for p in sources),encoding='utf-8')
            rows=['@echo off','call "'+str(setup)+'" intel64 >nul','if errorlevel 1 exit /b %errorlevel%',
                'cd /d "'+str(out)+'"','ifx /nologo /O2 /libs:dll /c "'+str(module)+'"','if errorlevel 1 exit /b %errorlevel%',
                'ifx /nologo /O2 /libs:dll /Qdiag-disable:8291 /c @sources.rsp','if errorlevel 1 exit /b %errorlevel%',
                'ifx /nologo /O2 /libs:dll /exe:"'+str(exe)+'" *.obj','exit /b %errorlevel%']
            (out/'build.bat').write_text('\r\n'.join(rows)+'\r\n',encoding='utf-8')
            result=subprocess.run(['cmd','/d','/c',str(out/'build.bat')],cwd=out,stdout=log,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW)
        else:
            command=shutil.which('gfortran')
            if not command:raise ValueError('Install gfortran first')
            flags=['-O2','-std=legacy','-fallow-argument-mismatch','-ffixed-line-length-none','-ffree-line-length-none']
            result=subprocess.run([command,*flags,'-c',str(module)],cwd=out,stdout=log,stderr=subprocess.STDOUT)
            if not result.returncode:
                short=[x for x in flags if x!='-ffixed-line-length-none']+['-ffixed-line-length-72']
                result=subprocess.run([command,*short,'-c',str(src/'tran.f')],cwd=out,stdout=log,stderr=subprocess.STDOUT)
            if not result.returncode:
                rest=[p for p in sources if p.name!='tran.f']
                result=subprocess.run([command,*flags,str(out/'modparm.o'),str(out/'tran.o'),*map(str,rest),'-o',str(exe)],
                    cwd=out,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode or not exe.is_file():
        print('\n'.join(log_path.read_text(errors='replace').splitlines()[-30:]),file=sys.stderr)
        raise RuntimeError('SWAT compilation failed; see '+str(log_path))
    print('[BUILD] '+str(exe),flush=True)
    return exe


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--population',type=int,default=POPULATION)
    p.add_argument('--generations',type=int,default=GENERATIONS)
    p.add_argument('--workers',type=int,default=WORKERS)
    p.add_argument('--seed',type=int,default=SEED)
    p.add_argument('--txtinout',type=Path,default=ROOT/'TxtInOut')
    p.add_argument('--exe',type=Path,default=ROOT/'build'/('SWAT_DPS.exe' if os.name=='nt' else 'SWAT_DPS'))
    p.add_argument('--output',type=Path,default=ROOT/'results'/'run')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--build',action='store_true')
    p.add_argument('--build-only',action='store_true')
    p.add_argument('--compiler',choices=['ifx','gfortran'])
    args=p.parse_args(argv)
    try:
        if args.build or args.build_only or not args.exe.exists():
            built=build_swat(args.compiler)
            if args.exe.resolve()!=(ROOT/'build'/('SWAT_DPS.exe' if os.name=='nt' else 'SWAT_DPS')).resolve():
                raise ValueError('Build writes ROOT/build/SWAT_DPS; omit custom --exe when building')
            args.exe=built
        if not args.build_only:run(args)
    except (ValueError,RuntimeError) as e:p.exit(1,'ERROR: '+str(e)+'\n')


if __name__=='__main__':main()
