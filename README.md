# Monte Carlo Simulation of the 2D Ising Model

A Fortran-based Monte Carlo simulator of the two-dimensional Ising model, used to study
the thermodynamic evolution of the spin lattice, compare update mechanisms (Metropolis,
Glauber, Wolff), analyse their ergodic / autocorrelation behaviour, and locate the
critical temperature of the ferromagnetic phase transition via finite-size scaling.

This is the codebase behind the paper *"Estudio Monte Carlo del modelo de Ising 2D"*
(`articulo/`, `articulo_en/` for the English version, compiled PDFs at
`article.pdf` / `artículo.pdf`) and an accompanying presentation (`presentacion/`).

## What the simulation does

The Ising model places a spin `s_i = ±1` on every site of an `L × L` square lattice
with periodic boundary conditions. Neighbouring spins interact ferromagnetically
(Hamiltonian `H = -J Σ s_i s_j`, nearest neighbours only, `J = 1`). At low temperature
the system prefers to align into large ordered domains (spontaneous magnetization); at
high temperature thermal noise destroys the order. In the infinite-lattice limit this
is a genuine continuous phase transition, with an exactly known critical temperature
(Onsager's solution): `T_c = 2.2691853 J/k_B`.

Because the partition function cannot be summed exactly for anything but small
lattices, the equilibrium averages (energy, magnetization, susceptibility, heat
capacity...) are estimated with Markov-chain Monte Carlo: the simulator builds a chain
of spin configurations whose stationary distribution is the Boltzmann distribution at
a given temperature, and time-averages over that chain once it has thermalized.

## Monte Carlo methods implemented

All three algorithms are implemented from scratch in a single Fortran file
(`simulator.f`, ~720 lines) and share the same lattice, temperature grid and
observable-accumulation code:

- **Metropolis**: pick a random spin, compute the energy cost `ΔE` of flipping it, and
  accept the flip with probability `min(1, e^{-ΔE/T})`. One Monte Carlo "step" here is
  `N = L²` such single-spin trial updates (a full sweep).
- **Glauber (heat-bath) dynamics**: pick a random spin and set it according to the
  local equilibrium probability `p = 1 / (1 + e^{ΔE/T})`, regardless of its current
  state. Same asymptotic (Boltzmann) distribution as Metropolis but different
  relaxation dynamics.
- **Wolff cluster algorithm**: instead of flipping one spin at a time, grow a cluster
  of aligned neighbouring spins (each bond added with probability `p_add = 1 - e^{-2J/T}`)
  and flip the whole cluster at once. This performs much better near `T_c`, where
  single-spin algorithms suffer from *critical slowing down* (very long
  autocorrelation times because large domains take a long time to flip spin-by-spin).

Both single-spin algorithms use a Metropolis-style local-field lookup table
(`exponentials(...)`, indexed by the discrete set of possible `ΔE` values on a square
lattice) instead of calling `exp()` inside the hot loop, and a hand-rolled
`xorshift64` PRNG instead of Fortran's built-in generator, for speed. The Wolff
routine grows clusters with an explicit stack (avoiding recursion) and only clears the
`visited` flags that were actually touched, rather than the whole lattice.

For every temperature the code records, once the chain has thermalized:

- mean energy per spin `⟨e⟩`
- mean absolute magnetization per spin `⟨|m|⟩`
- heat capacity per spin `c`, from energy fluctuations: `c = (⟨E²⟩-⟨E⟩²)/(N T²)`
- `⟨m²⟩`, `⟨m⁴⟩` (used downstream for the susceptibility and the Binder cumulant)
- the structure factor at the smallest lattice wavevector `⟨|m(k_min)|²⟩` (used to
  estimate the correlation length)
- the integrated autocorrelation time `τ_int` of `|m(t)|`, via the Madras–Sokal
  windowing method — this is the ergodic/mixing-time diagnostic mentioned on the CV:
  it directly measures how many Monte Carlo steps are needed between statistically
  independent samples, and how that scales with `L` and algorithm choice near `T_c`.

The temperature grid is not uniform: it is coarse away from `T_c` and refined close to
it (window width shrinking with `L`), so computational effort is concentrated where
the physics changes fastest.

Downstream (in `data_analysis.ipynb`), two more quantities are derived from the raw
columns:
- **Susceptibility** `χ = (L²/T)(⟨m²⟩ - ⟨|m|⟩²)`, which diverges at `T_c` as `L → ∞`.
- **Binder cumulant** `U₄ = 1 - ⟨m⁴⟩/(3⟨m²⟩²)`, whose curves for different `L` cross
  at (an estimate of) `T_c`, independent of finite-size corrections to leading order —
  the standard way to pin down a critical point from finite lattices.
- **Correlation length** `ξ_L = (1/(2 sin(π/L))) √(⟨m²⟩/⟨|m(k₁)|²⟩ - 1)`.

## What I actually ran in this session

I did **not** fabricate any numbers or images. Everything below is either (a) freshly
computed by me in this environment, or (b) the production data/figures that already
shipped in this repository, clearly labelled as such.

**Compiled and ran the simulator myself:**
```
gfortran -O3 -march=native -ffast-math simulator.f -o simulator
bash run.sh --smoke        # L = 16, 32; all three algorithms; 20000 MC steps/T
```
This completed in **2m 58s** on a 4-core Intel i7-10875H, and produced fresh,
independently generated data confirming the simulator builds and behaves physically:
magnetization drops from 1 to ~0 across the expected temperature window, energy rises
monotonically, and the drop sharpens as `L` grows from 16 to 32 — exactly the expected
finite-size trend of a continuous phase transition. That raw output is kept at
`results/smoke_test_raw_data/` and plotted in `results/smoke_test_verification.png`.

**Did not reproduce the full production run.** The repository already ships
high-statistics data for `L = 16, 32, 64, 128` (`datos/{metropolis,glauber,wolff}_2d.txt`,
~85-89 temperatures each, 20000 MC steps with 6000 discarded for thermalization). A
back-of-envelope timing extrapolation from my own `L=16` run (10s single-threaded)
puts a full `L=128` run at roughly 10-15 minutes *per algorithm* running alone, and the
repo's own `run.sh` launches all 12 (algorithm × L) jobs concurrently on what is only a
4-core machine here — so a full non-smoke run would run several jobs at a time per
core and plausibly take the better part of an hour. I judged that too slow/wasteful to
redo blindly in this session when the committed data is already the real output of
this exact simulator binary, so instead I:

1. Verified the binary that produces that data still compiles cleanly and runs
   correctly (the smoke test above).
2. **Recomputed every plot and every reported number below directly from the
   already-committed raw data in `datos/`**, using my own script
   (`results/generate_plots.py`, logic ported faithfully from the repo's own
   `data_analysis.ipynb`), rather than reusing the repo's pre-rendered SVGs in
   `Gráficas/`. The numbers below are therefore freshly computed by me, from real
   simulation output, not copied from the paper or notebook.

If you want the full `L=16,32,64,128` grid regenerated from scratch, run:
```
./run.sh            # full run — expect roughly 30-60+ minutes on a 4-core laptop
```

## Results

All figures below were generated in this session by `results/generate_plots.py`
running against `datos/*.txt` (except `smoke_test_verification.png`, generated from
the fresh smoke-test run described above).

### Energy and magnetization vs. temperature

![Energy vs Temperature](results/energy_vs_temperature.png)

![Magnetization vs Temperature](results/magnetization_vs_temperature.png)

Both observables show the expected finite-size behaviour: away from `T_c` all lattice
sizes agree closely; near `T_c` the transition sharpens visibly as `L` grows from 16 to
128, consistent with the true transition becoming a true discontinuity only in the
`L → ∞` limit. All three algorithms (Metropolis, Glauber, Wolff) agree with each other
within statistical noise, as they must — they sample the same equilibrium Boltzmann
distribution, only the dynamics used to get there differ.

### Heat capacity and susceptibility

![Heat capacity vs Temperature](results/heat_capacity_vs_temperature.png)

![Susceptibility vs Temperature](results/susceptibility_vs_temperature.png)

Both quantities peak near `T_c`, and the peak grows and sharpens with `L` — the
finite-size signature of the divergences (`c`, `χ → ∞`) predicted for the infinite
lattice.

### Critical temperature via Binder cumulant crossing

![Binder cumulant crossing](results/binder_cumulant_critical_temperature.png)

The Binder cumulant `U₄(T)` curves for different `L` cross close to `T_c`, largely
independent of `L`. Using pairwise crossings between consecutive lattice sizes and
extrapolating linearly in `1/L²` to the `L → ∞` limit (exactly as implemented in
`data_analysis.ipynb`, reproduced independently in `results/generate_plots.py`) gives,
**from the data in this repository**:

| Algorithm  | Extrapolated `T_c` (`L→∞`) | Mean of all pairwise crossings |
|------------|----------------------------|---------------------------------|
| Metropolis | 2.27259                    | 2.27007                          |
| Glauber    | 2.26932                    | 2.26995                          |
| Wolff      | 2.27196                    | 2.27437                          |
| **Onsager (exact)** | **2.26919**        | —                                |

All three estimates land within ~0.15% of the exact Onsager value, which is the kind
of agreement expected from `L up to 128` with the crossing method. The raw numbers are
saved at `results/critical_temperature_summary.json`.

### Simulator self-verification (this session)

![Smoke test verification](results/smoke_test_verification.png)

Fresh `L = 16, 32` run I compiled and executed myself in this session (not
pre-existing repo data), confirming the simulator reproduces the qualitative Ising
transition for all three algorithms before trusting the larger pre-existing dataset
above.

## Repository layout

- `simulator.f` — the Fortran simulator (Metropolis, Glauber, Wolff).
- `run.sh` — compiles and launches the full (or `--smoke`) parameter sweep.
- `datos/` — production simulation output (`L = 16,32,64,128`), as shipped in the repo.
- `data_analysis.ipynb` — the original analysis notebook (Tc, critical exponents,
  finite-size-scaling collapses).
- `results/` — plots and analysis generated in this session (see above), plus
  `generate_plots.py` (the script that produced them) and the raw smoke-test output.
- `Gráficas/`, `articulo/`, `articulo_en/`, `presentacion/` — the original paper,
  its English translation, the slide deck, and their pre-rendered figures (SVG/PDF),
  as shipped in the repo before this session.
- `article.pdf` / `artículo.pdf` — compiled versions of the paper.

## Requirements

- `gfortran` to compile `simulator.f` (verified with GNU Fortran 11.4.0 on Ubuntu 22.04).
- Python 3 with `numpy`, `pandas`, `matplotlib`, `scipy` (see `requirements.txt`) to
  run the analysis notebook / `results/generate_plots.py`.
- `latexmk` + a TeX distribution to rebuild the paper (not exercised in this session).
