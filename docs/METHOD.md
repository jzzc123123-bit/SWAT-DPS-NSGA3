# Annual fertilizer feedback and NSGA-III

SWAT-DPS-NSGA3 implements the eight-state fertilizer-management formulation
with **NSGA-III**. A single RBF policy is shared across subbasins.

## Beginning-of-year state

A single policy is shared across subbasins, but each subbasin supplies its own
agricultural-HRU-area-weighted state vector. Soil nitrogen pools are first summed
over soil layers within each HRU, then weighted by HRU area within the subbasin.
The eight inputs, in file order, are:

| Input | SWAT variable | Before normalization |
|---|---|---|
| Soil nitrate N | sum of `sol_no3` over layers | kg N/ha |
| Soil ammonium N | sum of `sol_nh3` over layers | kg N/ha |
| Fresh organic N | sum of `sol_fon` over layers | kg N/ha |
| Active organic N | sum of `sol_aorgn` over layers | kg N/ha |
| Stable organic N | sum of `sol_orgn` over layers | kg N/ha |
| Shallow-aquifer nitrate N | `shallst_n` | kg N/ha |
| Relative soil water | `sol_sw / sol_sumfc` | dimensionless |
| Shallow-aquifer water | `shallst` | mm |

The model calls `DPS_CAPTURE_YEAR_START` after `sim_inityr` and before the daily
loop, then calls `DPS_YEAR_START`. No future-year outcome is a policy input.
Bounds are fixed externally before optimization; out-of-range components are
logged and clipped only at the controller input, not in SWAT's physical stores.
The state vector is a reduced representation, not a proof of strict or sufficient
Markov information.

## Gaussian policy and action mapping

With normalized state `z=clip((x-lo)/(hi-lo),0,1)`, the raw policy is

```text
phi_i = exp(-sum_j ((z_j - c_ij)/b_ij)^2)
raw_k = a_k + sum_i w_ik * phi_i
u_k = clip(raw_k,0,1)
```

There is no factor1/2 in the exponent, no normalization by the sum of basis
activations, and no independent231-parameter policy per subbasin.

- `B=8`, `K=3`, `N=B+K+1=12`.
- `P=K+N(2B+K)=231` shared optimized parameters.
- Biases and centers: `[0,1]`; widths: `[0.01,1]`; weights: `[-1,1]`.
- Serialization order: biases, row-major weights `[N,K]`, centers `[N,B]`,
  widths `[N,B]`. This flat order is an interface choice, not a change in the RBF.

```text
annual_scale = clip(1.2 * sum(u), 0.2, 1.2)
trial_shares = u / sum(u)
shares = Euclidean projection of trial_shares onto:
         sum(shares)=1, lower=[0.3,0.1,0.3], upper=[0.5,0.3,0.5]
target_phase_N = fixed_baseline_annual_N * annual_scale * shares
phase_multiplier = target_phase_N / fixed_baseline_phase_N
```

For all-zero clipped signals, the allocation is40:20:40 and the annual multiplier
is0.2. The annual reference is the **fixed baseline**, not the previous year's
dose. The final shares are not independently clipped and renormalized: they are
projected jointly onto the bounded sum-one simplex.

Each agricultural HRU's existing fertilizer product amount is multiplied by
the corresponding subbasin/event multiplier. This retains within-subbasin
baseline differences. Only the three scheduled chemical-fertilizer operations
are controlled. Their dates, manure inputs and domestic inputs remain fixed.
In the published fork, urea is fertilizer database ID4. The example state bounds
are illustrative: users must establish and freeze appropriate bounds for their
own SWAT project rather than copying them into a scientific application.

## Four physical objectives and two strict constraints

Let `D[t,s,k]` be the executed phase N rate (kg N/ha), `D0[s]` the fixed baseline
annual N rate, `A[s]` agricultural area and `T` the number of evaluation years.
The minimized objective vector is:

1. Mean annual outlet TN, `sum(monthly_TN_kg)/T/1000`, t N/year.
2. Negative mean annual basin grain, `-sum(annual_subbasin_grain_t)/T`, dry t/year.
3. Action change in percentage points/year:
   `100/(T-1) * sum_t sum_s (A[s]/sum(A))*sum_k abs(D[t,s,k]-D[t-1,s,k])/D0[s]`.
4. Negative minimum subbasin cumulative grain retention:
   `-min_s sum_t Y[t,s]/sum_t Y0[t,s]`.

The smoothness objective uses **executed N rates**, not changes in raw RBF outputs.
It excludes warmup and the warmup-to-evaluation transition. For2008–2050 this
means43 evaluation years and42 transitions.

Feasibility requires:

- Long-term basin grain retention **strictly greater than0.95**.
- The fraction of months with outlet TN **strictly above4mg/L** is **strictly
  less than0.25**. For516 months this permits at most128 exceedances.

The constraint does not impose a maximum concentration for every month, an
annual grain floor, or a per-subbasin hard grain floor. The4mg/L threshold is a
research setting, not a statutory compliance assertion. Monthly concentration
is `TN_kg *1000 / monthly_outflow_volume_m3`; zero flow with positive TN fails
validation, while zero flow with zero TN is recorded as no outflow.

## Search and convergence reporting

The default full search is100 individuals,1400 generations including
initialization, and8 isolated workers. NSGA-III uses Das–Dennis reference
directions (84 directions at population100), SBX probability0.9/eta15, and
polynomial mutation probability1/per-variable probability`1/P`/eta20.

The four physical objectives are scaled by paired baseline TN and grain,100pp,
and unit local retention for the portable optimizer. Scaling and HV reference
are frozen before a run and recorded in its manifest. The minimal case uses
`[1.2,0,1.2,0]` in these minimized coordinates. This is not the numerical affine
scale/reference used to reconstruct the paper's historical HV. Historical HV
numbers therefore cannot be compared directly with this example. Set one
common fixed profile before comparing methods; never revise it mid-run.

`history.csv` reports feasible-archive HV and the number of policies outside the
fixed reference. Empty feasible archives legitimately have HV0. A fixed-budget
stop, a stable best objective or a nonzero HV does not prove convergence.
Checkpoints include both global RNG states and the algorithm state for correct
continuation with pymoo0.6.1.x and0.6.2.
