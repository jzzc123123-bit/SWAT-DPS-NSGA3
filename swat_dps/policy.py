"""Three-event nitrogen policy and independent Gaussian RBF/budget oracle."""
from pathlib import Path
import numpy as np
B,K,N,P=8,3,12,231
SCHEMA='DPS_SUBBASIN_SPLITN_V5'
STATES=('soil_nitrate_N','soil_ammonium_N','soil_fresh_org_N',
        'soil_active_org_N','soil_stable_org_N','aquifer_NO3_N',
        'soil_water','shallow_aquifer_water')
PHASES=('basal','topdress1','topdress2')
REGIONAL_SHARES=np.array([.4,.2,.4])
PHASE_LO=np.array([.3,.1,.3])
PHASE_HI=np.array([.5,.3,.5])
def vector(bias=None,weights=None,centers=None,widths=None):
    if bias is None:bias=REGIONAL_SHARES/1.2
    return np.r_[np.asarray(bias).ravel(),
        np.zeros((N,K)).ravel() if weights is None else np.asarray(weights).ravel(),
        np.full((N,B),.5).ravel() if centers is None else np.asarray(centers).ravel(),
        np.ones((N,B)).ravel() if widths is None else np.asarray(widths).ravel()]
def unpack(theta):
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
def evaluate(theta,state,lo,hi,baseline_phase_N,annual_bounds=(.2,1.2),phase_lower=PHASE_LO,phase_upper=PHASE_HI,phase_reference=REGIONAL_SHARES):
    a,w,c,r=unpack(theta)
    state=np.asarray(state);lo=np.asarray(lo);hi=np.asarray(hi)
    if state.shape!=lo.shape or lo.shape!=hi.shape or state.shape!=(B,):raise ValueError('Eight state scales required')
    if not np.isfinite(np.r_[state,lo,hi]).all() or not (hi>lo).all():raise ValueError('Invalid state scales')
    z=np.clip((state-lo)/(hi-lo),0,1)
    phi=np.exp(-np.sum(((z-c)/r)**2,axis=1))
    raw=a+phi@w;u=np.clip(raw,0,1)
    result=allocate(u,baseline_phase_N,annual_bounds,phase_lower,phase_upper,phase_reference)
    result.update(z=z,phi=phi,raw=raw,u=u)
    return result
def write(path,theta,lo,hi,mode=1,start=2008,annual_bounds=(.2,1.2),phase_lower=PHASE_LO,phase_upper=PHASE_HI,phase_reference=REGIONAL_SHARES,outlet=21):
    a,w,c,r=unpack(theta)
    evaluate(theta,lo,lo,hi,np.ones(3),annual_bounds,phase_lower,phase_upper,phase_reference)
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
def read(path):
    rows=Path(path).read_text(encoding='ascii').splitlines()
    if len(rows)!=46 or rows[0]!=SCHEMA:raise ValueError('Expected SPLITN V5, 46 rows')
    if tuple(map(int,rows[2].split()))!=(B,K,N,P):raise ValueError('Policy dimensions')
    mode,start,outlet=map(int,rows[1].split())
    expected=[2,3,3,3,8,8,3]+[3]*N+[8]*(2*N)
    if any(len(row.split())!=width for row,width in zip(rows[3:],expected)):raise ValueError('Policy row width')
    bounds,phase_lower,phase_upper,phase_reference,lo,hi=[np.array(list(map(float,row.split()))) for row in rows[3:9]]
    theta=np.array([float(x) for row in rows[9:] for x in row.split()]);unpack(theta)
    evaluate(theta,lo,lo,hi,np.ones(3),bounds,phase_lower,phase_upper,phase_reference)
    if mode not in (0,1) or outlet<1:raise ValueError('Policy mode/outlet')
    return dict(mode=mode,start=start,outlet=outlet,bounds=bounds,lo=lo,hi=hi,theta=theta,
                phase_lower=phase_lower,phase_upper=phase_upper,phase_reference=phase_reference)
