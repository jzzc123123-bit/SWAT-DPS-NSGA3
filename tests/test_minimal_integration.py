import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path
import numpy as np
from swat_dps import policy
from swat_dps.minimal import MinimalCase
from swat_dps.metrics import evaluate
from swat_dps.optimize import load_config

ROOT=Path(__file__).resolve().parents[1]

class IntegrationTests(unittest.TestCase):
    def test_local_states_change_actions_and_policy_changes_future_states(self):
        c=load_config(ROOT/'examples/minimal_case/config.json');model=MinimalCase(c)
        bias=[.25,.15,.25];w=np.zeros((12,3));w[0]=[.08,-.03,.06]
        centers=np.full((12,8),.3);widths=np.full((12,8),.8)
        theta=policy.vector(bias,w,centers,widths)
        t=model.simulate(theta);base=model.simulate(theta,baseline=True)
        self.assertEqual(t.phase_N.shape,(6,3,3));self.assertEqual(t.monthly_TN_kg.size,72)
        self.assertFalse(np.allclose(t.phase_N[:,0]/t.baseline_phase_N[0],t.phase_N[:,1]/t.baseline_phase_N[1]))
        self.assertFalse(np.allclose(t.phase_N[0]/t.baseline_phase_N,t.phase_N[-1]/t.baseline_phase_N))
        self.assertFalse(np.allclose(t.monthly_TN_kg,base.monthly_TN_kg))
        self.assertEqual(evaluate(t,base,c)['transitions'],5)

    def launch(self,path,generations,resume=False):
        command=[sys.executable,str(ROOT/'swat_nsga3.py'),'--config',str(ROOT/'examples/minimal_case/config.json'),
                 '--population','8','--generations',str(generations),'--workers','2','--seed','12','--output',str(path)]
        if resume:command.append('--resume')
        done=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,timeout=120)
        self.assertEqual(done.returncode,0,done.stdout+'\n'+done.stderr)
        return json.loads((path/'status.json').read_text())

    def test_resume_matches_uninterrupted_search_and_exports_feasible_policy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=root/'resumed';b=root/'fresh'
            self.launch(a,2);status=self.launch(a,4,True);fresh=self.launch(b,4)
            self.assertEqual(status['evaluation_requests'],32);self.assertGreater(status['feasible_archive_size'],0)
            self.assertAlmostEqual(status['HV'],fresh['HV'],places=12)
            with np.load(a/'pareto_policies.npz') as aa,np.load(b/'pareto_policies.npz') as bb:
                self.assertEqual(set(aa['keys']),set(bb['keys']))
            self.assertTrue((a/'example_policy/actions.csv').is_file())
            self.assertTrue(json.loads((a/'example_policy_metrics.json').read_text())['feasible'])
            with (a/'history.csv').open() as stream:rows=list(csv.DictReader(stream))
            self.assertEqual([int(r['generation']) for r in rows],[1,2,3,4])
            self.assertTrue(all(float(rows[i]['HV'])<=float(rows[i+1]['HV'])+1e-12 for i in range(3)))

if __name__=='__main__':unittest.main()
