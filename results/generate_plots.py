"""
Reproduces the core analysis of data_analysis.ipynb against the repo's
production Monte Carlo data (datos/*.txt), and generates a smaller set of
independent PNG plots + a critical-temperature estimate, saved under
results/. Logic (Tc via Binder-cumulant crossing, chi, xi definitions) is
copied faithfully from the notebook so the numbers are directly comparable.

This script does not fabricate any numbers: it recomputes everything from
the raw simulation output already committed in datos/.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from itertools import combinations

plt.rcParams.update({
    'figure.figsize':    (7, 5),
    'font.size':         12,
    'axes.titlesize':    13,
    'axes.labelsize':    12,
    'xtick.labelsize':   10,
    'ytick.labelsize':   10,
    'legend.fontsize':   9,
    'legend.framealpha':  0.85,
    'legend.edgecolor':  '0.8',
    'lines.linewidth':   1.4,
    'lines.markersize':  3,
    'axes.grid':         True,
    'grid.alpha':        0.3,
    'grid.linestyle':    '--',
    'axes.spines.top':   False,
    'axes.spines.right': False,
})

TC_EXACT = 2.2691853  # Onsager exact solution

ALGOS = {
    'Metropolis': 'datos/metropolis_2d.txt',
    'Glauber':    'datos/glauber_2d.txt',
    'Wolff':      'datos/wolff_2d.txt',
}
COLS = ['L', 'T', 'e', 'm', 'c', 'm2', 'm4', 'mk2', 'tau']

L_COLORS  = {16: '#e41a1c', 32: '#377eb8', 64: '#4daf4a', 128: '#984ea3'}
L_MARKERS = {16: 'o',       32: 's',       64: '^',       128: 'D'}


def _load(fname):
    df = pd.read_csv(fname, sep=r'\s+', header=None, names=COLS)
    m2_conn = df['m2'] - df['m']**2
    df['U4']  = 1.0 - df['m4'] / (3.0 * df['m2']**2)
    df['chi'] = (df['L']**2 / df['T']) * m2_conn
    df['xi']  = (np.sqrt(np.maximum(m2_conn / df['mk2'] - 1.0, 0.0))
                 / (2.0 * np.sin(np.pi / df['L'])))
    return df


def slice_L(df, L):
    return df[df['L'] == L].sort_values('T')


def find_Tc(df, L_values, t_window=(1.8, 2.8)):
    interps = {L: PchipInterpolator(slice_L(df, L)['T'].values,
                                     slice_L(df, L)['U4'].values,
                                     extrapolate=False)
               for L in L_values}

    def _find_crossings(pairs):
        crossings, inv_L2 = [], []
        for L1, L2 in pairs:
            s1, s2 = slice_L(df, L1), slice_L(df, L2)
            T_lo = max(s1['T'].min(), s2['T'].min(), t_window[0])
            T_hi = min(s1['T'].max(), s2['T'].max(), t_window[1])
            Ts   = np.linspace(T_lo, T_hi, 2000)
            diff = interps[L1](Ts) - interps[L2](Ts)
            roots = []
            for j in range(len(diff) - 1):
                if diff[j] == 0:
                    roots.append(Ts[j])
                elif diff[j] * diff[j+1] < 0:
                    roots.append(brentq(
                        lambda t, _j=j: interps[L1](t) - interps[L2](t),
                        Ts[j], Ts[j+1],
                    ))
            if roots:
                crossings.append(float(np.median(roots)))
                inv_L2.append(1.0 / L2**2)
        return np.array(crossings), np.array(inv_L2)

    all_crossings, _ = _find_crossings(combinations(L_values, 2))
    tc_mean = float(np.mean(all_crossings)) if len(all_crossings) else np.nan

    consec_crossings, inv_L2 = _find_crossings(zip(L_values[:-1], L_values[1:]))
    if len(consec_crossings) >= 2:
        slope, tc_inf = np.polyfit(inv_L2, consec_crossings, 1)
    elif len(consec_crossings) == 1:
        slope, tc_inf = np.nan, float(consec_crossings[0])
    else:
        slope, tc_inf = np.nan, np.nan

    return {
        'Tc_mean': tc_mean,
        'Tc_inf': float(tc_inf),
        'slope': float(slope) if not np.isnan(slope) else np.nan,
        'crossings': consec_crossings,
        'inv_L2': inv_L2,
    }


data = {algo: _load(fname) for algo, fname in ALGOS.items()}
L_values = sorted(data['Glauber']['L'].unique())
print("L values found:", L_values)

# ---------------------------------------------------------------
# 1) Energy per spin vs T
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
for ax, (algo, df) in zip(axes, data.items()):
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['e'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=3, lw=1.2, label=f'L={L}')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.2, label='Onsager $T_c$')
    ax.set_xlabel('T'); ax.set_title(algo); ax.set_xlim(0, 5)
axes[0].set_ylabel(r'Energy per spin $\langle e \rangle$')
axes[0].legend(fontsize=8)
fig.suptitle('Ising 2D — Energy per spin vs Temperature (Metropolis / Glauber / Wolff)')
fig.tight_layout()
fig.savefig('results/energy_vs_temperature.png', dpi=150)
plt.close(fig)

# ---------------------------------------------------------------
# 2) Magnetization per spin vs T
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
for ax, (algo, df) in zip(axes, data.items()):
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['m'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=3, lw=1.2, label=f'L={L}')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.2, label='Onsager $T_c$')
    ax.set_xlabel('T'); ax.set_title(algo); ax.set_xlim(0, 5)
axes[0].set_ylabel(r'Magnetization per spin $\langle |m| \rangle$')
axes[0].legend(fontsize=8)
fig.suptitle('Ising 2D — Magnetization per spin vs Temperature')
fig.tight_layout()
fig.savefig('results/magnetization_vs_temperature.png', dpi=150)
plt.close(fig)

# ---------------------------------------------------------------
# 3) Heat capacity vs T
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
for ax, (algo, df) in zip(axes, data.items()):
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['c'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=3, lw=1.2, label=f'L={L}')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.2, label='Onsager $T_c$')
    ax.set_xlabel('T'); ax.set_title(algo); ax.set_xlim(1.5, 3.5)
axes[0].set_ylabel(r'Heat capacity per spin $c$')
axes[0].legend(fontsize=8)
fig.suptitle('Ising 2D — Heat capacity vs Temperature (peak sharpens near $T_c$ as L grows)')
fig.tight_layout()
fig.savefig('results/heat_capacity_vs_temperature.png', dpi=150)
plt.close(fig)

# ---------------------------------------------------------------
# 4) Susceptibility vs T
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
for ax, (algo, df) in zip(axes, data.items()):
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['chi'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=3, lw=1.2, label=f'L={L}')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.2, label='Onsager $T_c$')
    ax.set_xlabel('T'); ax.set_title(algo); ax.set_xlim(1.5, 3.5)
axes[0].set_ylabel(r'Susceptibility $\chi$')
axes[0].legend(fontsize=8)
fig.suptitle(r'Ising 2D — Magnetic susceptibility vs Temperature (diverges at $T_c$ as $L\to\infty$)')
fig.tight_layout()
fig.savefig('results/susceptibility_vs_temperature.png', dpi=150)
plt.close(fig)

# ---------------------------------------------------------------
# 5) Binder cumulant crossing -> Tc
# ---------------------------------------------------------------
fits = {algo: find_Tc(df, L_values) for algo, df in data.items()}

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
T_ZOOM = (2.15, 2.40)
for ax, (algo, df) in zip(axes, data.items()):
    for L in L_values:
        s = slice_L(df, L)
        s_zoom = s[(s['T'] >= T_ZOOM[0]) & (s['T'] <= T_ZOOM[1])]
        ax.plot(s_zoom['T'], s_zoom['U4'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=3, lw=1.4, label=f'L={L}')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.2, label='Onsager $T_c$')
    Tc = fits[algo]['Tc_inf']
    if not np.isnan(Tc):
        ax.axvline(Tc, color='k', ls='-', lw=1.0, alpha=0.7,
                   label=f'Estimated $T_c$={Tc:.4f}')
    ax.set_xlabel('T'); ax.set_title(algo)
axes[0].set_ylabel(r'Binder cumulant $U_4$')
axes[0].legend(fontsize=8)
fig.suptitle(r'Binder cumulant $U_4(T)$ crossing $\Rightarrow$ critical temperature estimate')
fig.tight_layout()
fig.savefig('results/binder_cumulant_critical_temperature.png', dpi=150)
plt.close(fig)

# ---------------------------------------------------------------
# Print Tc summary (used verbatim in README / report — not fabricated)
# ---------------------------------------------------------------
summary = {}
for algo, df in data.items():
    res = fits[algo]
    Tc  = res['Tc_inf'] if not np.isnan(res['Tc_inf']) else res['Tc_mean']
    summary[algo] = {
        'Tc_extrapolated_Linfty': round(Tc, 5),
        'Tc_mean_all_pairs': round(res['Tc_mean'], 5),
    }
summary['Onsager_exact'] = TC_EXACT
print(json.dumps(summary, indent=2))
with open('results/critical_temperature_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("Plots written to results/")
