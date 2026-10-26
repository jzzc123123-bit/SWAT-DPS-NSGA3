# SWAT-DPS-NSGA3

Complete modified classic-SWAT source, a runnable minimum TxtInOut case,
and one Python script for eight-state RBF-DPS optimization with NSGA-III.

```text
SWAT_Source/       complete modified SWAT2016 Rev.664 source (312 files)
TxtInOut/          minimum example: 2 subbasins, 2 agricultural HRUs
NSGA_III.py        compilation, SWAT execution and NSGA-III optimization
requirements.txt  Python dependencies
```

## Run

Use Python 3.10-3.12. Install Intel oneAPI Fortran on Windows, or gfortran
on Linux. Run these commands from the repository root:

```sh
python -m pip install -r requirements.txt
python NSGA_III.py --build-only
python NSGA_III.py --population 12 --generations 3 --workers 2
```

Every evaluation invokes the compiled SWAT executable. There is no surrogate
or simplified replacement model. Compilation writes to `build/`; optimization
writes to `results/run/`. These generated files are not distributed.

For a larger search, use a new output directory:

```sh
python NSGA_III.py --population 100 --generations 1400 --workers 8 --output results/full
```

Add `--resume` with the same settings/output to continue a checkpoint.
`Ctrl+C` requests a stop after the current generation. Console progress and
`history.csv` report feasibility and HV. `pareto.csv`/`pareto_policies.npz`
contain feasible nondominated solutions; `example_policy/` contains a replay
of the minimum-TN feasible policy found. This is not a paper-policy selection.

## Settings and method

Edit the constants and `CASE` dictionary near the top of `NSGA_III.py`.
Population, generations, workers and seed also have command-line options.
For your own project, change `CASE` years/outlet, eight fixed state ranges
and 27 fixed hydrological values; pass `--txtinout` and `--exe` as needed.
The supplied values are illustrative, not calibrated study-basin parameters.

The shared policy has B=8 states, K=3 outputs, N=12 RBFs and 231 parameters.
Beginning-of-year agricultural-HRU states are area-weighted per subbasin.
Annual chemical N is 20%-120% of its fixed baseline. Basal/topdress1/topdress2
shares are jointly constrained to 30%-50%/10%-30%/30%-50%, summing to one.
Manure and scheduled application dates remain fixed.

Objectives: minimize annual outlet TN, maximize annual basin grain,
minimize interannual changes in executed phase N rates, and maximize the
minimum subbasin cumulative grain retention. Constraints: basin grain >95%
of baseline; outlet TN >4 mg/L in <25% of evaluation months.
The example HV reference is `[1.2, 0, 1.2, 0]` in baseline-relative coordinates.

## Minimum case

The case runs 2000-2005, with 2000-2001 warmup and 2002-2005 evaluation.
It has scheduled maize planting/harvest, three urea applications per year,
fixed manure and phosphorus inputs, and independently specified daily weather.
It contains no study-basin observations, fitted policies or optimization results.
It demonstrates execution, not calibration or reproduction of paper results.

Input layouts and standard database excerpts follow the public
[SWATdata demo](https://github.com/chrisschuerz/SWATdata).
Attribution and the applicable third-party terms are in `TxtInOut/NOTICE.txt`
and `TxtInOut/COPYING.txt`. This does not assign a new license to the SWAT-DPS code.

The source's parameter-audit loop uses the actual subbasin count; physical
process equations and the eight-state DPS contract are retained.
