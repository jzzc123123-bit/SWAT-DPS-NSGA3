"""Build the current SWAT fork or a small oracle using the actual Fortran routines."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',choices=['ifx','gfortran'],default='ifx' if os.name=='nt' else 'gfortran')
    parser.add_argument('--oracle',action='store_true')
    parser.add_argument('--output',type=Path,default=Path('build'))
    args=parser.parse_args();root=Path(__file__).resolve().parents[1];src=root/'swat_source/dps_b8k3'
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    suffix='.exe' if os.name=='nt' else '';name='policy_oracle' if args.oracle else 'SWAT_DPS_B8K3'
    exe=out/(name+suffix)
    if args.oracle:sources=[src/'dps_rbf.f',src/'dps_feedback.f',root/'tools/policy_oracle.f90']
    else:sources=sorted(p for p in src.iterdir() if p.suffix.lower() in {'.f','.f90'} and p.name!='modparm.f')
    module=src/'modparm.f'
    if args.compiler=='ifx':
        if os.name!='nt':raise SystemExit('This helper uses the Windows ifx environment; choose gfortran on Linux')
        oneapi=Path(os.environ.get('ONEAPI_ROOT',r'C:\Program Files (x86)\Intel\oneAPI'))/'setvars.bat'
        if not oneapi.is_file():raise SystemExit('Set ONEAPI_ROOT to your Intel oneAPI installation')
        rsp=out/'sources.rsp';rsp.write_text('\n'.join('"'+str(p)+'"' for p in sources),encoding='utf-8')
        batch=out/'build.bat'
        rows=['@echo off','setlocal','call "'+str(oneapi)+'" intel64 >nul','if errorlevel 1 exit /b %errorlevel%',
              'cd /d "'+str(out)+'"','ifx /nologo /O2 /libs:dll /c "'+str(module)+'"','if errorlevel 1 exit /b %errorlevel%',
              'ifx /nologo /O2 /libs:dll /Qdiag-disable:8291 /c @sources.rsp','if errorlevel 1 exit /b %errorlevel%',
              'ifx /nologo /O2 /libs:dll /exe:"'+str(exe)+'" *.obj','exit /b %errorlevel%']
        batch.write_text('\r\n'.join(rows)+'\r\n',encoding='utf-8')
        cmd=['cmd.exe','/d','/c',str(batch)]
        with (out/'compile.log').open('w',encoding='utf-8') as log:
            done=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
    else:
        compiler=shutil.which('gfortran')
        if not compiler:raise SystemExit('gfortran is not installed or not on PATH')
        flags=['-O2','-std=legacy','-fallow-argument-mismatch','-ffixed-line-length-none','-ffree-line-length-none']
        # tran.f alone places a redundant ampersand beyond column72. Its
        # executable statements fit within72 columns; compile that unchanged
        # native file at its intended width rather than editing its equations.
        short_form=[p for p in sources if p.name=='tran.f']
        rest=[p for p in sources if p.name!='tran.f']
        short_objects=[]
        with (out/'compile.log').open('w',encoding='utf-8') as log:
            done=subprocess.run([compiler,*flags,'-c',str(module)],cwd=out,stdout=log,stderr=subprocess.STDOUT)
            for p in short_form:
                if done.returncode:break
                short_flags=[x for x in flags if x!='-ffixed-line-length-none']+['-ffixed-line-length-72']
                done=subprocess.run([compiler,*short_flags,'-c',str(p)],cwd=out,stdout=log,stderr=subprocess.STDOUT)
                short_objects.append(str(out/(p.stem+'.o')))
            if not done.returncode:done=subprocess.run([compiler,*flags,str(out/'modparm.o'),*short_objects,*map(str,rest),'-o',str(exe)],cwd=out,stdout=log,stderr=subprocess.STDOUT)
    if done.returncode or not exe.is_file():
        lines=(out/'compile.log').read_text(encoding='utf-8',errors='replace').splitlines()
        print('\n'.join(lines[-100:]),file=sys.stderr)
        raise SystemExit('Build failed; see '+str(out/'compile.log'))
    manifest={'compiler':args.compiler,'oracle':args.oracle,'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),
              'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [module,*sources]}}
    (out/'build_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(exe)

if __name__=='__main__':main()
