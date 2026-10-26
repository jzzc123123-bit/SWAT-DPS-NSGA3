"""Small executable case for the policy/optimizer interface, not SWAT physics.

All forcing and initial conditions are defined here or in the example JSON.
No watershed observations, calibrated project or paper results are used.
"""
import numpy as np
from . import policy
from .metrics import Trajectory

class MinimalCase:
    name='minimal'
    def __init__(self,config,workspace=None):
        self.config=config;self.c=config['minimal_case']
        self.membership=np.asarray(self.c['subbasins'],dtype=int)
        self.area=np.asarray(self.c['HRU_area_ha'],dtype=float)
        self.base=np.asarray(self.c['baseline_N_kg_ha'],dtype=float)
        if self.area.shape!=self.membership.shape or self.base.shape!=self.area.shape or (self.area<=0).any() or (self.base<=0).any():
            raise ValueError('Minimal case HRU geometry')
        self.subs=np.unique(self.membership);self.subarea=np.array([self.area[self.membership==s].sum() for s in self.subs])
        self.shares=np.asarray(self.c['baseline_phase_shares'],dtype=float)
        self.base_phase=np.array([np.average(self.base[self.membership==s,None]*self.shares,axis=0,weights=self.area[self.membership==s]) for s in self.subs])

    def fingerprint(self):
        return {'backend':'minimal','configuration':self.config}

    def simulate(self,theta,baseline=False):
        c=self.config;first,last=map(int,c['evaluation_years']);warm_first,warm_last=map(int,c['warmup_years'])
        if warm_last!=first-1:raise ValueError('Warmup must end immediately before evaluation')
        years=np.arange(first,last+1);ny=len(years);ns=len(self.subs);nh=len(self.area)
        # Soil N: NO3, NH4, fresh, active, stable; shallow aquifer nitrate,
        # relative soil water and shallow aquifer water. Area-weighted at Jan1.
        h=np.arange(nh,dtype=float)
        state=np.c_[30+2*h,np.full(nh,1.),70+3*h,240+5*h,2600+20*h,12+h,np.full(nh,.7),150+5*h]
        phase=np.zeros((ny,ns,3));grain=np.zeros((ny,ns));loads=np.zeros((ny,12));volumes=np.zeros((ny,12));oob=0
        lo=np.asarray(c['state_lower']);hi=np.asarray(c['state_upper'])
        rainfall_shape=np.array([.025,.025,.045,.07,.12,.16,.17,.15,.11,.06,.035,.03])
        rainfall_shape/=rainfall_shape.sum()
        uptake_shape=np.array([0,0,0,.06,.15,.23,.25,.2,.11,0,0,0])
        for y in range(warm_first,last+1):
            doses=self.base[:,None]*self.shares
            if y>=first and not baseline:
                for j,s in enumerate(self.subs):
                    mask=self.membership==s;x=np.average(state[mask],axis=0,weights=self.area[mask])
                    oob+=int(np.count_nonzero((x<lo)|(x>hi)))
                    out=policy.evaluate(theta,x,lo,hi,self.base_phase[j],c['annual_bounds'],c['phase_lower'],c['phase_upper'],c['phase_reference'])
                    doses[mask]*=out['phase_scale']
            if y>=first:
                iy=y-first
                for j,s in enumerate(self.subs):
                    mask=self.membership==s;phase[iy,j]=np.average(doses[mask],axis=0,weights=self.area[mask])
            taken=np.zeros(nh);water_score=np.zeros(nh)
            for month in range(1,13):
                precip=self.c['annual_precipitation_mm']*rainfall_shape[month-1]*(1+.12*np.sin(.9*(y-warm_first)))
                precip=precip*(1+.03*h)
                bucket=state[:,6]*180+precip
                evap=28+15*np.sin(np.pi*(month-1)/11)
                drain=np.maximum(0,bucket-180);bucket=np.clip(bucket-evap-drain,5,180)
                state[:,6]=bucket/180
                for k,m in enumerate(self.c['application_months']):
                    if month==m:state[:,1]+=doses[:,k] # mineral N as NH4, urea
                if month==3:state[:,2]+=self.c['fixed_manure_N_kg_ha_year']
                fresh=state[:,2]*.04;state[:,2]-=fresh;state[:,3]+=fresh*.5;state[:,1]+=fresh*.5
                active=state[:,3]*.004;state[:,3]-=active;state[:,1]+=active
                stable=state[:,4]*.0002;state[:,4]-=stable;state[:,3]+=stable
                nitrified=state[:,1]*.38;state[:,1]-=nitrified;state[:,0]+=nitrified
                available=state[:,0]+state[:,1]
                uptake=np.minimum(available,180*uptake_shape[month-1]*(.65+.35*state[:,6]))
                fraction=np.divide(state[:,0],available,out=np.zeros_like(available),where=available>0)
                state[:,0]-=uptake*fraction;state[:,1]-=uptake*(1-fraction);taken+=uptake
                leached=state[:,0]*np.minimum(.35,drain/350);state[:,0]-=leached;state[:,5]+=leached
                runoff=state[:,0]*np.minimum(.03,precip/12000);state[:,0]-=runoff
                denitrified=state[:,0]*.006;state[:,0]-=denitrified
                gw_release=state[:,5]*.065;state[:,5]-=gw_release
                state[:,7]=np.maximum(5,state[:,7]+drain*.6-state[:,7]*.04)
                local_export=(runoff+gw_release)*self.area
                # A three-reach chain; upstream export experiences more attenuation.
                routing=np.array([.9**(ns-int(np.where(self.subs==s)[0][0])-1) for s in self.membership])
                mass=float(np.sum(local_export*routing)+self.c['fixed_domestic_N_kg_month'])
                volume=float(np.sum((drain*.6+state[:,7]*.04)*self.area*10)+self.c['background_flow_volume_m3_month'])
                water_score+=state[:,6]/12
                if y>=first:loads[iy,month-1]=mass;volumes[iy,month-1]=max(volume,1.)
            yield_t_ha=8*(1-np.exp(-taken/60))*(.8+.2*water_score)
            state[:,2]+=12+2*yield_t_ha # residues remain in the soil for next year
            if y>=first:
                for j,s in enumerate(self.subs):
                    mask=self.membership==s;grain[iy,j]=np.sum(yield_t_ha[mask]*self.area[mask])
        return Trajectory(years,self.subs,self.subarea,self.base_phase,phase,grain,loads,volumes,oob).validate()
