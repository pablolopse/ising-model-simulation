"""Generate individual per-algorithm figures for the Beamer presentation."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
os.chdir(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.interpolate import interp1d
from shutil import copy2

OUTDIR = 'presentacion/figs'

plt.rcParams.update({
    'font.size': 11, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 8,
    'legend.framealpha': 0.85, 'legend.edgecolor': '0.8',
    'lines.linewidth': 1.4, 'lines.markersize': 2,
    'axes.grid': True, 'grid.alpha': 0.3, 'grid.linestyle': '--',
    'axes.spines.top': False, 'axes.spines.right': False,
})

TC_EXACT = 2.2691853
ALGOS = {'Metropolis': 'datos/metropolis_2d.txt',
         'Glauber':    'datos/glauber_2d.txt',
         'Wolff':      'datos/wolff_2d.txt'}
COLS = ['L', 'T', 'e', 'm', 'c', 'm2', 'm4', 'mk2', 'tau']
L_COLORS  = {16: '#e41a1c', 32: '#377eb8', 64: '#4daf4a', 128: '#984ea3'}
L_MARKERS = {16: 'o', 32: 's', 64: '^', 128: 'D'}

def _load(fname):
    df = pd.read_csv(fname, sep=r'\s+', header=None, names=COLS)
    m2_conn = df['m2'] - df['m']**2
    df['U4']  = 1.0 - df['m4'] / (3.0 * df['m2']**2)
    df['chi'] = (df['L']**2 / df['T']) * m2_conn
    df['xi']  = (np.sqrt(np.maximum(m2_conn / df['mk2'] - 1.0, 0.0))
                 / (2.0 * np.sin(np.pi / df['L'])))
    return df

data = {algo: _load(fname) for algo, fname in ALGOS.items()}
L_values = sorted(data['Glauber']['L'].unique())

def slice_L(df, L):
    return df[df['L'] == L].sort_values('T')

SZ = (4.5, 3.4)  # single-panel size — fits in 1/3 of a 16:9 slide

def savefig(name):
    path = f'{OUTDIR}/{name}'
    plt.savefig(path, bbox_inches='tight')
    plt.close()
    print(f'  {path}')

# ── c ──────────────────────────────────────────────────────────────────────
for algo, df in data.items():
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['c'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=2, lw=1.4, label=f'$L={L}$')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.0, label=r'$T_c$')
    ax.set_xlim(1.5, 3.5)
    ax.set_xlabel(r'$T$'); ax.set_ylabel(r'$c$')
    ax.set_title(algo); ax.legend(loc='upper right', fontsize=7)
    plt.tight_layout()
    savefig(f'c_{algo}.svg')

# ── m ──────────────────────────────────────────────────────────────────────
for algo, df in data.items():
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['m'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=2, lw=1.4, label=f'$L={L}$')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.0, label=r'$T_c$')
    ax.set_xlim(1.5, 3.5)
    ax.set_xlabel(r'$T$'); ax.set_ylabel(r'$\langle|m|\rangle$')
    ax.set_title(algo); ax.legend(loc='upper right', fontsize=7)
    plt.tight_layout()
    savefig(f'm_{algo}.svg')

# ── fss_binder ─────────────────────────────────────────────────────────────
T_ZOOM = (2.20, 2.35)
for algo, df in data.items():
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        s_zoom = s[(s['T'] >= T_ZOOM[0]) & (s['T'] <= T_ZOOM[1])]
        ax.plot(s_zoom['T'], s_zoom['U4'], color=L_COLORS[L],
                marker=L_MARKERS[L], ms=2, lw=1.4, label=f'$L={L}$')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.0, label=r'$T_c$')
    ax.set_xlabel(r'$T$'); ax.set_ylabel(r'$U_4$')
    ax.set_title(algo); ax.legend(loc='lower left', fontsize=7)
    plt.tight_layout()
    savefig(f'fss_binder_{algo}.svg')

# ── suscept ────────────────────────────────────────────────────────────────
T_CHI = (1.9, 2.7)
for algo, df in data.items():
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['chi'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=2, lw=1.4, label=f'$L={L}$')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.0, label=r'$T_c$')
    ax.set_xlim(*T_CHI)
    ax.set_xlabel(r'$T$'); ax.set_ylabel(r'$\chi$')
    ax.set_title(algo); ax.legend(fontsize=7)
    plt.tight_layout()
    savefig(f'suscept_{algo}.svg')

# ── suscept_fss ─────────────────────────────────────────────────────────────
FSS_PARAMS = {
    'Metropolis': {'Tc': 2.273, 'nu': 0.975, 'gamma': 1.7},
    'Glauber':    {'Tc': 2.269, 'nu': 0.995, 'gamma': 1.76},
    'Wolff':      {'Tc': 2.272, 'nu': 0.985, 'gamma': 1.74},
}
for algo, params in FSS_PARAMS.items():
    Tc, nu, gamma = params['Tc'], params['nu'], params['gamma']
    df = data[algo]
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        x = ((s['T'].values - Tc) / Tc) * (L ** (1.0 / nu))
        y = s['chi'].values * (L ** (-gamma / nu))
        ax.plot(x, y, color=L_COLORS[L], marker=L_MARKERS[L],
                ms=2, alpha=0.8, lw=1.2, label=f'$L={L}$')
    ax.set_title(f'{algo}\n'
                 f'$T_c={Tc:.3f}$, $\\nu={nu:.3f}$, $\\gamma={gamma:.3f}$',
                 fontsize=9)
    ax.set_ylabel(r'$\chi\,L^{-\gamma/\nu}$')
    ax.set_xlabel(r'$(T-T_c)/T_c\cdot L^{1/\nu}$')
    ax.set_xlim(-1, 3); ax.legend(fontsize=7)
    plt.tight_layout()
    savefig(f'suscept_fss_{algo}.svg')

# ── xi_vs_T ────────────────────────────────────────────────────────────────
for algo, df in data.items():
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['xi'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=2, lw=1.4, label=f'$L={L}$')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.0, label=r'$T_c$')
    ax.set_xlabel(r'$T$'); ax.set_ylabel(r'$\xi_L$')
    ax.set_title(algo); ax.legend(fontsize=7)
    plt.tight_layout()
    savefig(f'xi_vs_T_{algo}.svg')

# ── xi_ratio ────────────────────────────────────────────────────────────────
for algo, df in data.items():
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['xi'] / L, color=L_COLORS[L], marker=L_MARKERS[L],
                ms=2, lw=1.4, label=f'$L={L}$')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.0, label=r'$T_c$')
    ax.set_xlabel(r'$T$'); ax.set_ylabel(r'$\xi_L / L$')
    ax.set_title(algo); ax.legend(fontsize=7)
    plt.tight_layout()
    savefig(f'xi_ratio_{algo}.svg')

# ── xi_fss ──────────────────────────────────────────────────────────────────
for algo, params in FSS_PARAMS.items():
    Tc, nu = params['Tc'], params['nu']
    df = data[algo]
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        x = ((s['T'].values - Tc) / Tc) * (L ** (1.0 / nu))
        y = s['xi'].values / L
        ax.plot(x, y, color=L_COLORS[L], marker=L_MARKERS[L],
                ms=2, alpha=0.8, lw=1.2, label=f'$L={L}$')
    ax.set_title(f'{algo}\n$T_c={Tc:.3f}$, $\\nu={nu:.3f}$', fontsize=9)
    ax.set_ylabel(r'$\xi_L / L$')
    ax.set_xlabel(r'$(T-T_c)/T_c\cdot L^{1/\nu}$')
    ax.set_xlim(-1, 3); ax.set_ylim(bottom=0); ax.legend(fontsize=7)
    plt.tight_layout()
    savefig(f'xi_fss_{algo}.svg')

# ── tau_vs_T ────────────────────────────────────────────────────────────────
for algo, df in data.items():
    fig, ax = plt.subplots(figsize=SZ)
    for L in L_values:
        s = slice_L(df, L)
        ax.plot(s['T'], s['tau'], color=L_COLORS[L], marker=L_MARKERS[L],
                ms=2, lw=1.4, label=f'$L={L}$')
    ax.axvline(TC_EXACT, color='gray', ls=':', lw=1.0, label=r'$T_c$')
    ax.set_xlabel(r'$T$'); ax.set_ylabel(r'$\tau_\mathrm{int}$')
    ax.set_title(algo); ax.legend(fontsize=7)
    plt.tight_layout()
    savefig(f'tau_vs_T_{algo}.svg')

# ── chi_scaling, xi_scaling, tau_scaling, extrp_lin (already 1×3 from article)
for name in ['chi_scaling', 'xi_scaling', 'tau_scaling', 'extrp_lin']:
    copy2(f'article/figs/{name}.svg', f'{OUTDIR}/{name}.svg')
    print(f'  {OUTDIR}/{name}.svg (copied)')

print('\nDone.')
