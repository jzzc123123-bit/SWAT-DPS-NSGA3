"""Independent Python active-set projection versus the actual Fortran policy."""
from pathlib import Path
import argparse,subprocess,sys,tempfile,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from swat_dps import policy
from swat_dps.swat_backend import number

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--cases',type=int,default=300);args=parser.parse_args()
    rng=np.random.default_rng(41026);lo=np.zeros(8);hi=np.array([500,10,300,600,6000,80,1.5,600]);rows=[];expected=[]
    for i in range(args.cases):
        a=rng.random(3);w=rng.uniform(-1,1,(12,3));centers=rng.random((12,8));widths=rng.uniform(.01,1,(12,8))
        if i<7:
            a=np.array([[0,0,0],[1,1,1],[1,0,0],[0,1,0],[0,0,1],[.01,.02,.03],[.2,.1,.3]][i],dtype=float);w[:]=0
        theta=policy.vector(a,w,centers,widths);state=rng.uniform(-.1,1.3,8)*hi;base=np.exp(rng.uniform(1,6,3))
        expected.append(policy.evaluate(theta,state,lo,hi,base))
        for x in [lo,hi,state,base,a,*w,*centers,*widths]:rows.append(' '.join(format(float(v),'.17g') for v in x))
    with tempfile.TemporaryDirectory(prefix='swat_dps_oracle_') as temp:
        p=Path(temp);(p/'oracle_input.dat').write_text(str(args.cases)+'\n'+'\n'.join(rows)+'\n',encoding='ascii')
        subprocess.run([str(args.executable.resolve())],cwd=p,check=True,capture_output=True)
        values=np.loadtxt(p/'dps_actions.tsv',converters=number)
        if values.shape!=(args.cases,65):raise ValueError('Oracle output dimensions')
        for i,out in enumerate(expected):
            np.testing.assert_allclose(values[i,12],out['annual_scale'],atol=1e-12)
            for sl,key in [(slice(19,22),'target'),(slice(28,31),'shares'),(slice(39,47),'z'),(slice(47,50),'raw'),(slice(50,53),'u'),(slice(53,65),'phi')]:
                np.testing.assert_allclose(values[i,sl],out[key],atol=1e-10,rtol=1e-10,err_msg=str((i,key)))
    print(json.dumps({'passed':True,'cases':args.cases,'checked':'actual Fortran RBF and annual action mapping; Python independent simplex projection'}))

if __name__=='__main__':main()
