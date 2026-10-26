import copy,json,math,tempfile,unittest
from pathlib import Path
import numpy as np
from swat_dps import policy
from swat_dps.metrics import Trajectory,evaluate
from swat_dps.optimize import bounds,load_config
from swat_dps.swat_backend import number,SwatCase

ROOT=Path(__file__).resolve().parents[1]

class MethodTests(unittest.TestCase):
    def setUp(self):self.config=load_config(ROOT/'examples/minimal_case/config.json')

    def test_dimensions_and_parameter_bounds(self):
        self.assertEqual(policy.N,policy.B+policy.K+1)
        self.assertEqual(policy.P,policy.K+policy.N*(2*policy.B+policy.K))
        lo,hi=bounds();self.assertEqual(len(lo),231)
        self.assertTrue(np.all(lo[:3]==0));self.assertTrue(np.all(lo[3:39]==-1));self.assertTrue(np.all(lo[-96:]==.01))
        self.assertTrue(np.all(hi==1))

    def test_zero_signals_preserve_annual_floor_and_reference_shares(self):
        result=policy.allocate(np.zeros(3),np.array([80,60,20]))
        self.assertAlmostEqual(result['annual_scale'],.2)
        np.testing.assert_allclose(result['target'],32*np.array([.4,.2,.4]))

    def test_projection_preserves_mass_and_bounds(self):
        # Hand-solvable Euclidean projection of [1,0,0] onto these bounds.
        out=policy.allocate([1,0,0],[80,60,20])
        np.testing.assert_allclose(out['shares'],[.5,.2,.3],atol=1e-12)
        self.assertAlmostEqual(sum(out['target']),1.2*160)
        # Spatial baseline differences survive the shared event multiplier.
        baseline=np.array([[80,60,20],[160,120,40]])
        realized=baseline*out['phase_scale']
        np.testing.assert_allclose(realized[1],2*realized[0])

    def test_policy_roundtrip_supports_a_different_outlet(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'dps_policy.in';v=policy.vector()
            policy.write(p,v,np.zeros(8),np.ones(8),outlet=3)
            result=policy.read(p);self.assertEqual(result['outlet'],3)
            np.testing.assert_array_equal(result['theta'],v)

    def fixture(self):
        years=np.array([2002,2003]);subs=np.array([1,2]);area=np.array([1.,3.]);base=np.tile([30.,10.,60.],(2,1))
        baseline=Trajectory(years,subs,area,base,np.tile(base,(2,1,1)),np.tile([10.,30.],(2,1)),np.full((2,12),2.),np.full((2,12),1000.))
        current=Trajectory(years,subs,area,base,np.array([base,[[20.,20.,50.],[40.,10.,60.]]]),np.array([[9.,31.],[9.,29.]]),np.ones((2,12)),np.full((2,12),1000.))
        return current,baseline

    def test_objectives_against_hand_calculation(self):
        current,base=self.fixture();out=evaluate(current,base,self.config)
        np.testing.assert_allclose(out['F'],[.012,-39.,15.,-.9])
        self.assertAlmostEqual(out['grain_retention'],.975);self.assertEqual(out['transitions'],1)
        self.assertEqual(out['worst_subbasin'],1);self.assertTrue(out['feasible'])

    def test_strict_thresholds_and_concentration_units(self):
        current,base=self.fixture()
        current.grain_t=base.grain_t*.95
        self.assertGreater(evaluate(current,base,self.config)['G'][0],0)
        current.grain_t=base.grain_t.copy()
        current.monthly_TN_kg[:]=4.0 # 4kg N in1000m3 =4mg/L; equality not exceedance
        self.assertEqual(evaluate(current,base,self.config)['TN_exceedance_count'],0)
        current.monthly_TN_kg.flat[:5]=4.1
        self.assertTrue(evaluate(current,base,self.config)['feasible'])
        current.monthly_TN_kg.flat[5]=4.1
        self.assertFalse(evaluate(current,base,self.config)['feasible'])
        self.assertEqual(math.ceil(.25*516)-1,128)

    def test_invalid_outputs_fail_instead_of_silent_repair(self):
        current,base=self.fixture();current.monthly_volume_m3[0,0]=0
        with self.assertRaises(ValueError):evaluate(current,base,self.config)
        for token in ['NaN','Infinity','********','not-a-number']:
            with self.assertRaises(ValueError):number(token)
        self.assertAlmostEqual(number('1.0864749155846127-171'),1.0864749155846127e-171)
        self.assertEqual(number(b'1.2D+01'),12.)

    def test_generated_native_dat_files_are_not_treated_as_inputs(self):
        for name in ['959999watout.dat','959999hru.dat','959999input.std','959999parameter_audit.txt','959999output.rch']:
            self.assertFalse(SwatCase._is_input(Path(name)),name)
        for name in ['1p.dat','21p.dat','000010001.hru','fert.dat','file.cio']:
            self.assertTrue(SwatCase._is_input(Path(name)),name)

if __name__=='__main__':unittest.main()
