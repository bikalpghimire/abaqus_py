# -*- coding: utf-8 -*-
"""Compare runs pairwise (coarse vs fine mesh), load solid / unload dashed, up to 35 mm.
Writes force_disp_<a>_vs_<b>.png for each pair."""
import glob
import os
import matplotlib.pyplot as plt
from force_disp import read_csv

HERE = os.path.dirname(os.path.abspath(__file__))
DMAX = 35.0
PAIRS = [
    ('All flexible', [('15', 'mesh 5 mm / 0.2 mm', 'tab:red'),
                      ('17', 'mesh 3 mm / 0.1 mm', 'tab:blue')]),
    ('Beam + connection flexible', [('14', 'mesh 5 mm / 0.2 mm', 'tab:red'),
                                    ('16', 'mesh 3 mm / 0.1 mm', 'tab:blue')]),
]

for title, runs in PAIRS:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for prefix, label, color in runs:
        data = read_csv(glob.glob(os.path.join(HERE, prefix + '-*_fd.csv'))[0])
        for step, ls in (('Load', '-'), ('Unload', '--')):
            if step not in data:
                continue
            u, f = data[step]
            keep = [(a, b) for a, b in zip(u, f) if a <= DMAX]
            ax.plot([k[0] for k in keep], [k[1] for k in keep], color=color, ls=ls, lw=1.3,
                    label='%s - %s' % (prefix, label) if step == 'Load' else None)
    ax.plot([], [], color='0.3', lw=1.3, label='loading')
    ax.plot([], [], color='0.3', lw=1.3, ls='--', label='unloading')
    ax.axhline(0.0, color='0.6', lw=0.6)
    ax.set_xlim(0, DMAX)
    ax.set_xlabel('Displacement of top plate [mm]')
    ax.set_ylabel('Reaction force [N]')
    ax.set_title('%s: run %s vs %s (up to %g mm)' % (title, runs[0][0], runs[1][0], DMAX))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(HERE, 'force_disp_%s_vs_%s.png' % (runs[0][0], runs[1][0]))
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print('Wrote ' + out)
