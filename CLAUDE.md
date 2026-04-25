# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TFG (Trabajo de Fin de Grado) computational physics project: Finite Size Scaling (FSS) analysis of the 2D Ising model using three Monte Carlo algorithms (Metropolis, Glauber, Wolff). The goal is to determine the critical temperature and critical exponents near the phase transition.

## Build & Run

**Compile and run the unified simulator** (generates all three `.txt` data files):
```bash
gfortran simulator.f -o simulator
./simulator
```

**Run analysis:**
```bash
jupyter notebook fss_analysis.ipynb
```

The standalone `*_2d.f` files and their binaries have been removed; `simulator.f` is the current single source of truth.

## Architecture

### Simulation layer (`simulator.f`)

Fortran 77-style program that runs all three algorithms sequentially for L = 16, 32, 64, 128.

**Output columns** (`metropolis_2d.txt`, `glauber_2d.txt`, `wolff_2d.txt`): `L  T  e  m  Cv  m2  m4`
- `e` = energy/spin, `m` = ⟨|magnetization|⟩/spin, `Cv` = heat capacity, `m2` = ⟨m²⟩, `m4` = ⟨m⁴⟩

**Temperature grid** (`BUILD_TEMPS`): adaptive per-L grid descending from T=5.0.
- Coarse step Δt=0.1 outside the fine window
- Fine step Δt=0.005 within `[Tc - 2/L, Tc + 3/L]` (the FSS-relevant range for scaled temperature X ∈ [−2, 3])

**MC step counts** (same rule for all three algorithms):
- Near Tc (T ∈ [2.0, 2.5]): STEPS=2,000,000, NSKIP=600,000
- Elsewhere: STEPS=20,000, NSKIP=6,000
- Metropolis/Glauber sample every 50 steps; Wolff samples every 5 sweeps (a sweep = enough cluster flips to turn ~N spins total)

**Key subroutines:**
- `INITIALIZE` — random ±1 lattice, computes initial E and M
- `BUILD_TEMPS` — constructs the adaptive temperature array
- `RUN_METROPOLIS` / `RUN_GLAUBER` / `RUN_WOLFF` — run one full temperature sweep and write output
- `METROPOLIS_MC_STEP` / `GLAUBER_MC_STEP` / `WOLFF_MC_STEP` — single MC step; precompute acceptance exponentials before the temperature loop to avoid repeated `exp()` calls

### Analysis layer (`fss_analysis.ipynb`)

Loads the three `.txt` files, adds derived columns `U4 = 1 - m4/(3·m2²)` (Binder cumulant) and `chi = L²·(m2 - m²)` (susceptibility), then:

1. **`find_Tc`** — estimates T_c from U₄ crossings between consecutive L pairs using linear interpolation + `brentq`; extrapolates to L→∞ via linear fit of T_cross vs 1/L².
2. **Polynomial crossing cell** — same logic using quadratic fits near T_c.
3. **`get_exponents`** — log-log fits to extract ν (from |dU₄/dT| at T_c), β/ν (from ⟨|m|⟩ at T_c), γ/ν (from χ_max).
4. **Interactive collapse widgets** (`ipywidgets`) — sliders for T_c, ν, γ to visually verify FSS collapse of χ for each algorithm.
5. **Extrapolation plots** — `Gráficas/extrp_lin.png`, `extrp_poli.png`; FSS plots saved to `Gráficas/fss_*.png`.

### Key physics parameters

- Coupling constant J = 1 (ferromagnetic), periodic boundary conditions
- Theoretical T_c = 2.2691853 (Onsager exact solution); search window for crossings: (1.8, 2.8)
- Exact exponents: ν = 1, β/ν = 1/8, γ/ν = 7/4
