# Using the complete SWAT backend

The minimal case is immediately executable without a watershed project. For a
real SWAT run, supply a complete independently prepared classic-SWAT TxtInOut
project, its calibrated parameter vector, and a build of the source in
`swat_source/dps_b8k3`. No research-basin input project or calibrated vector is
distributed in this release.

## Build

Windows, Intel oneAPI ifx (run from the repository root):

```powershell
python tools/build_swat.py --compiler ifx --output build/swat
```

The script discovers `ONEAPI_ROOT` or the standard oneAPI installation. The
current Windows build uses DLL runtime libraries to avoid static CRT symbol
collisions on recent Visual Studio versions. If runtime DLLs are not on PATH,
run from an initialized oneAPI command prompt or install the corresponding
Intel redistributables. Compiler logs and source/executable hashes go to the
ignored build directory. Do not mix oracle/full-build object files: choose
different output directories.

Linux/GNU build helper:

```sh
python tools/build_swat.py --compiler gfortran --output build/swat
```

This helper supplies legacy fixed-format and argument-mismatch flags. A
successful build on one compiler does not establish compatibility with every
toolchain. The repository workflow checks the full GNU build and policy oracle.

## Configure your project

Copy `examples/swat_case.template.json` outside version control and provide:

- Your `template_dir` and compiled `executable`.
- Exactly27 calibrated hydrological command-line values, in the existing
  `main.f` external order; these values are fixed during policy search.
- Fixed eight-dimensional state bounds appropriate to your own project.
- Simulation/warmup/evaluation years and the outlet reach number.

The illustrative state bounds in the template are **not** calibrated limits.
The chemical-fertilizer annual and phase bounds, constraints and four objective
definitions implement the current manuscript formulation. The template uses
the study's year layout only to illustrate file configuration, not as data.

The inherited main program takes27 calibration values and an output ID. The
external order is CANMX, relative CN2, CH_N2, CH_K2, relative SOL_AWC, relative
SOL_K, relative SOL_BD, CH_W2, CH_L2, ALPHA_BF, ESCO, GW_DELAY, CH_D, TLAPS,
TIMP, GWQMN, SMFMX, GW_REVAP, REVAPMN, CH_S2, SMFMN, SNO50COV, SNOCOVMX,
SOL_ALB, SURLAG, SFTMP, EPCO. Read the injection code and calibrate your own
model; this wrapper does not recalibrate these values. The original source's
calibration range guards remain present.

The backend verifies `file.cio` has monthly output (`IPRINT=0`), complete
calendar years and the matching `NYSKIP`. It expects exactly three positive
scheduled urea (ID4) operations and one annual `HARV/KILL` maize harvest per
agricultural HRU. It fails clearly for incomplete harvest/action/monthly
coverage or nonfinite output instead of silently filling missing values.
Projects with rotations or multiple harvests need an explicitly revised
harvest contract before using this reader.

```sh
python swat_nsga3.py --config /path/to/my_case.json \
  --population 100 --generations 1400 --workers 8 --output results/my_case
```

Every worker gets a separate copied TxtInOut; the original template is not
modified. Executable/input hashes, policy schema, environment and objective
settings are checked when resuming. Native monthly TN masses and flow volumes
determine the water-quality constraint. Native maize harvest kg/ha are
aggregated using independently parsed HRU areas to produce dry tonnes.

On POSIX systems the backend links the inherited Windows output filename
`NUL` to `/dev/null`, avoiding an unintended large HRU-output file. This changes
only where suppressed output is written, not the model's physical equations.

## Source changes and call chain

- `modparm.f`: shared dimensions, parameters, states and annual event scales.
- `simulate.f`: year-start capture and annual decisions before daily operations.
- `dps_feedback.f`: file parsing, agricultural mask, state aggregation, annual
  N budget and bounded-simplex allocation, trace outputs and outlet logging.
- `dps_rbf.f`: Gaussian basis evaluation and clipped control signals.
- `sched_mgt.f`: multiply the scheduled urea product amount by its event scale
  immediately before the existing fertilizer routine executes.
- `fert.f`: retains its existing fertilizer-to-soil-pool process calculations.
- `writea.f`: annual outlet diagnostic capture.

`dps_policy.in` is a46-line text input, schema`DPS_SUBBASIN_SPLITN_V5`; Python
and Fortran share the exact ordering. Mode0 retains baseline management and
mode1 activates annual RBF feedback from the configured start year. Warmup
management operates normally but is excluded from objectives. Missing policy
input disables feedback, but the optimizer always writes an explicit file.
The control is not implemented by changing point-source nitrogen inputs.

## Independent arithmetic cross-check

```sh
python tools/build_swat.py --oracle --output build/oracle
python tools/check_fortran_policy.py --executable build/oracle/policy_oracle.exe
```

On Linux omit`.exe`. The oracle calls the **actual published Fortran routines**
`DPS_RBF_EVAL` and`DPS_YEAR_START`. It checks random policies, clipping,
zero-output fallback and allocation against Python's independent active-set
simplex projection. This is an interface/arithmetic check, not a hydrological
calibration or a substitute for running a complete SWAT project.
