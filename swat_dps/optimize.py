"""Four-objective NSGA-III with isolated evaluations, fixed HV and checkpoints."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,csv,hashlib,json,math,os,platform,random,shutil,signal,sys,threading,time
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
from . import policy
from .metrics import evaluate,Trajectory,write_trajectory

def atomic_json(path,value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp.%s.%s'%(os.getpid(),threading.get_ident()))
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8');os.replace(tmp,p)

def key_for(theta):return hashlib.sha256(np.asarray(theta,dtype='<f8').tobytes()).hexdigest()

def bounds():
    return (np.r_[np.zeros(policy.K),-np.ones(policy.N*policy.K),np.zeros(policy.N*policy.B),np.full(policy.N*policy.B,.01)],
            np.ones(policy.P))

def initial_population(size,seed):
    rng=np.random.default_rng(seed);lo,hi=bounds();x=rng.uniform(lo,hi,(size,policy.P))
    # Generic constant policies, not optimized vectors from a research basin.
    for i,scale in enumerate([1.,.8,.6,.4,1.2,.2]):
        if i>=size:break
        x[i]=policy.vector(policy.REGIONAL_SHARES*scale/1.2)
    return x

def load_config(path):
    p=Path(path).resolve();c=json.loads(p.read_text(encoding='utf-8-sig'))
    for name in ['state_lower','state_upper']:
        if len(c[name])!=policy.B or not np.isfinite(c[name]).all():raise ValueError('Eight finite '+name+' values required')
    if not (np.asarray(c['state_upper'])>c['state_lower']).all():raise ValueError('Fixed state bounds must have positive widths')
    policy.evaluate(policy.vector(),c['state_lower'],c['state_lower'],c['state_upper'],np.ones(3),
                    c['annual_bounds'],c['phase_lower'],c['phase_upper'],c['phase_reference'])
    if list(c['annual_bounds'])!=[.2,1.2] or not np.allclose(c['phase_lower'],[.3,.1,.3]) or not np.allclose(c['phase_upper'],[.5,.3,.5]):
        raise ValueError('This release uses the manuscript annual and phase bounds')
    if not np.allclose(c['phase_reference'],[.4,.2,.4]):raise ValueError('Expected40:20:40 fallback')
    if (c.get('grain_floor_fraction',.95),c.get('TN_limit_mg_L',4.),c.get('TN_max_exceedance_fraction',.25))!=(.95,4.,.25):
        raise ValueError('This manuscript profile requires grain>95% and TN>4mg/L months<25%')
    if c['backend'] not in ['minimal','swat']:raise ValueError('Backend must be minimal or swat')
    if c['backend']=='swat':
        for field in ['template_dir','executable']:
            v=Path(c['swat'][field]).expanduser();c['swat'][field]=str(v if v.is_absolute() else (p.parent/v).resolve())
    ref=np.asarray(c['HV_reference_scaled'],dtype=float)
    if ref.shape!=(4,) or not np.isfinite(ref).all():raise ValueError('Four finite fixed HV reference values required')
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
    np.savez_compressed(output/'pareto_policies.npz',theta=np.asarray([r['theta'] for r in archive]).reshape(-1,policy.P),
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
    config=load_config(args.config);out=Path(args.output).resolve()
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
    if config['backend']=='minimal':
        from .minimal import MinimalCase
        backend=MinimalCase(config,runtime)
    else:
        from .swat_backend import SwatCase
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
        baseline=backend.simulate(policy.vector(),baseline=True);save_trajectory(baseline_path,baseline)
        write_trajectory(out/'baseline',baseline)
    baseline_hash=hashlib.sha256(baseline_path.read_bytes()).hexdigest()
    baseline_result=evaluate(baseline,baseline,config);atomic_json(out/'baseline_metrics.json',baseline_result)
    partitions=max(h for h in range(1,args.population) if math.comb(h+3,3)<=args.population)
    refs=get_reference_directions('das-dennis',4,n_partitions=partitions)
    xl,xu=bounds();problem=Problem(n_var=policy.P,n_obj=4,n_ieq_constr=2,xl=xl,xu=xu)
    if previous is None:
        algorithm=NSGA3(ref_dirs=refs,pop_size=args.population,sampling=initial_population(args.population,args.seed),
            crossover=SBX(prob=.9,eta=15),mutation=PM(prob=1.,prob_var=1/policy.P,eta=20),eliminate_duplicates=True)
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
        trajectory=backend.simulate(theta);result=evaluate(trajectory,baseline,config)
        result.update(key=key,signature=signature);atomic_json(p,result)
        return result,False
    manifest={'schema':'SWAT_DPS_NSGAIII_V2','signature':signature,'backend':config['backend'],'population':args.population,
        'generations':args.generations,'seed':args.seed,'workers':args.workers,'B':policy.B,'K':policy.K,'N':policy.N,'P':policy.P,
        'objectives':4,'constraints':2,'reference_directions':len(refs),'reference_partitions':partitions,
        'crossover':{'type':'SBX','probability':.9,'eta':15},'mutation':{'type':'PM','probability':1.,'per_variable':1/policy.P,'eta':20},
        'configuration':config,'runtime':{k:signature_data[k] for k in ['python','numpy','pymoo','dill']},
        'HV_reference_scaled':ref.tolist(),'HV_scope':'feasible cumulative nondominated archive, unchanged reference and scales',
        'generation_count_includes_initialization':True,'minimum_case_is_not_SWAT':config['backend']=='minimal'}
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
            policy.write(out/'dps_policy.in',theta,config['state_lower'],config['state_upper'],start=config['evaluation_years'][0],
                         annual_bounds=config['annual_bounds'],phase_lower=config['phase_lower'],phase_upper=config['phase_upper'],
                         phase_reference=config['phase_reference'],outlet=int(config['outlet_reach']))
            executed=backend.simulate(theta);result=evaluate(executed,baseline,config)
            np.testing.assert_allclose(result['F'],best['F'],rtol=1e-10,atol=1e-10)
            atomic_json(out/'example_policy_metrics.json',result);write_trajectory(out/'example_policy',executed)
        print('[%s] %s'%(status.upper(),out),flush=True)
    except Exception as e:
        atomic_json(out/'status.json',{'status':'failed','completed_generations':completed,'error':str(e)})
        raise
    finally:signal.signal(signal.SIGINT,old_signal)

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[1]/'examples/minimal_case/config.json')
    p.add_argument('--population',type=int,default=100);p.add_argument('--generations',type=int,default=1400)
    p.add_argument('--workers',type=int,default=8);p.add_argument('--seed',type=int,default=20261007)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--resume',action='store_true')
    args=p.parse_args(argv)
    try:run(args)
    except (ValueError,RuntimeError) as e:p.exit(1,'ERROR: '+str(e)+'\n')

if __name__=='__main__':main()
