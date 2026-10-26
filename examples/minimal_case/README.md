# Minimal example

From the repository root:

```sh
pip install -r requirements.txt
python examples/minimal_case/run.py
```

This runs20 individuals,5 generations and2 workers, and produces
`results/minimal_case`. It normally takes seconds. To change the budget:

```sh
python examples/minimal_case/run.py --population 100 --generations 30 \
  --workers 8 --output results/minimal_larger
```

The small case contains three subbasins and six agricultural HRUs, two warmup
years and six evaluation years. All inputs are defined by `config.json` and
the equations in `swat_dps/minimal.py`. It does not read a research watershed,
observations, calibrated parameters, real weather records, policy results or
files under the legacy data directory.

It exercises the same231-parameter policy, local HRU-area aggregation, annual
N/phase mapping, four objectives, two strict constraints, NSGA-III, feasible
archive, fixed-reference HV and checkpoint interfaces. Its small monthly
soil/water/N/crop process model is an executable interface example, **not
the SWAT process model**, a calibrated case or a reproduction of paper results.
The complete SWAT source and backend are provided separately and require the
user's own input project. See `docs/SWAT_INTEGRATION.md`.

Main outputs:

- `history.csv`: per-generation HV, feasibility, evaluation/cache counts and time.
- `pareto.csv` and`pareto_policies.npz`: all feasible nondominated solutions found.
- `example_policy.json` and`dps_policy.in`: the minimum-TN feasible example
  policy, not the paper's selected representative policy.
- `example_policy/actions.csv`, `annual_grain.csv`, `monthly_outlet.csv`:
  re-evaluated example-policy outputs.
- `baseline/`: paired baseline outputs from the same small process model.
- `manifest.json`, `status.json`, `checkpoint.pkl`: run settings and recovery.

Resume without restarting completed evaluations:

```sh
python examples/minimal_case/run.py --generations 10 --resume
```

Keep the same population, seed, configuration, code and environment. A target
generation count may increase; workers may change. Only load trusted local
checkpoints. `Ctrl+C` requests a stop after the current generation; a file named
`STOP` in the output directory requests the same. Remove`STOP` before resuming.
Only one optimizer may own an output directory. If the process was forcibly
terminated, check the PID in`optimizer.lock` before removing a stale lock.
