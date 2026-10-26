"""The four manuscript objectives, evaluated on executed actions and paired outputs."""
from dataclasses import dataclass
import math
import numpy as np

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

def evaluate(current,baseline,config):
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
