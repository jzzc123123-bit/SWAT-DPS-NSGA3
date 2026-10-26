# SWAT-DPS-NSGA3

Annual fertilizer feedback using a shared eight-state Gaussian RBF policy
embedded in classic SWAT2016 Rev.664, optimized with **NSGA-III**.

This repository contains the current implementation and an independently
specified **minimal example**. It does not distribute research-basin inputs,
observations, calibrated parameter vectors, fitted policies or study results.

## Run the minimal example

From the repository root, using Python 3.10-3.12:

```sh
python -m pip install -r requirements.txt
python examples/minimal_case/run.py
```

The default run uses 20 individuals, 5 generations and 2 workers. It writes
results to `results/minimal_case` and normally completes in seconds. The case
contains 3 subbasins, 6 agricultural HRUs, 2 warmup years and 6 evaluation years.

The example uses a small monthly soil/water/nitrogen/crop model to demonstrate
the policy, action mapping, objectives, constraints and optimizer end to end.
**It does not invoke SWAT or reproduce the paper's numerical results.** The
complete SWAT source and execution backend are supplied separately for users
with their own complete classic-SWAT input project.

Change the search budget or use a new output directory:

```sh
python examples/minimal_case/run.py --population 100 --generations 100 --workers 8 --output results/my_example
```

See [minimal-example instructions](examples/minimal_case/README.md).

## Policy and optimization

| Component | Implementation |
|---|---|
| Decision time | Beginning of each year, before daily simulation |
| State aggregation | Agricultural-HRU area weighting within each subbasin |
| Eight states | Soil NO3, NH4, fresh/active/stable organic N; shallow-aquifer NO3; relative soil water; shallow-aquifer water |
| Shared RBF policy | B=8, K=3, N=B+K+1=12, P=K+N(2B+K)=231 |
| Spatial action | One shared rule; each subbasin's states determine its actions |
| Annual chemical N | 20%-120% of the fixed baseline annual N rate |
| Phase shares | Basal 30%-50%, topdress1 10%-30%, topdress2 30%-50%; sum=1 |
| Fixed management | Chemical application dates, manure and domestic inputs |
| Full-search defaults | NSGA-III, population 100, 1400 generations, 8 workers |

Four objectives: minimize mean annual outlet TN load; maximize mean annual
basin grain production; minimize area-weighted interannual changes in executed
phase N rates; maximize the minimum subbasin long-term grain retention.

Two strict constraints: long-term basin grain retention >95%, and fewer than
25% of months with outlet TN >4 mg/L. The concentration threshold is a research
setting. See [method and formulas](docs/METHOD.md).

## Run a user-supplied SWAT project

Build the source with Intel oneAPI ifx on Windows or gfortran on Linux:

```sh
python tools/build_swat.py --output build/swat
```

Copy `examples/swat_case.template.json` to your own configuration file and
provide your complete TxtInOut path, executable, 27 fixed hydrological values,
outlet and fixed state bounds. The template's bounds are illustrative.

```sh
python swat_nsga3.py --config /path/to/my_case.json --population 100 --generations 1400 --workers 8 --output results/my_swat_run
```

Workers use separate input copies. See [SWAT build and integration](docs/SWAT_INTEGRATION.md)
for required fertilizer operations, native outputs and the source call chain.

## Monitoring and recovery

The terminal reports generation, evaluations, feasible policies, HV and elapsed
time. Outputs include `history.csv`, `pareto.csv`, `pareto_policies.npz`,
`manifest.json`, `status.json`, `checkpoint.pkl`, `dps_policy.in`, and baseline
and re-evaluated policy tables. Outputs and build artifacts are ignored by Git.

HV uses fixed baseline-relative objective coordinates and a fixed reference
recorded in the run manifest. The minimal example's reference is
`[1.2, 0, 1.2, 0]`; its values are not the paper's historical HV values.

`Ctrl+C` or a file named `STOP` requests a stop after the current generation.
Remove `STOP` and add `--resume` to the same command to continue. Keep the
population, seed, configuration, code and environment fixed; the target
generation count may increase. Only load trusted local checkpoints.

## Verification

```sh
python tools/check_release.py
python -m unittest discover -s tests -v
python tools/build_swat.py --oracle --output build/oracle
python tools/check_fortran_policy.py --executable build/oracle/policy_oracle.exe
```

On Linux omit `.exe`. The oracle calls the actual Fortran policy and annual
allocation routines. GitHub Actions runs Python tests, the minimal case,
Fortran/Python arithmetic checks and a full GNU SWAT build.

## Source layout

```text
swat_nsga3.py                    optimizer entry
swat_dps/                       policy, objectives, optimizer and model backends
swat_source/dps_b8k3/            complete current classic-SWAT source
examples/minimal_case/          runnable minimal example and configuration
examples/swat_case.template.json user-project configuration template
tools/                          build helpers and arithmetic/release checks
tests/                          method and optimizer integration tests
docs/                           method, integration and source-hash manifest
```
