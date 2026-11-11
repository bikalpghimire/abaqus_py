# -*- coding: utf-8 -*-
"""Force-displacement of runs 16, 17, 18 (load solid, unload dashed), clipped to 35 mm.
Writes one figure per run: force_disp_<run>.png"""
import glob
import os
import matplotlib.pyplot as plt
from force_disp import read_csv

HERE = os.path.dirname(os.path.abspath(__file__))
DMAX = 35.0
RUNS = [('16', 'Beam + connection flexible', 'tab:blue'),
        ('17', 'All flexible', 'tab:orange'),
        ('18', 'Beam flexible', 'tab:green')]

for prefix, label, color in RUNS:
    csv = glob.glob(os.path.join(HERE, prefix + '-*_fd.csv'))[0]
    data = read_csv(csv)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for step, ls in (('Load', '-'), ('Unload', '--')):
        if step not in data:
            continue
        u, f = data[step]
        keep = [(a, b) for a, b in zip(u, f) if a <= DMAX]
        ax.plot([k[0] for k in keep], [k[1] for k in keep], color=color, ls=ls, lw=1.3,
                label='loading' if step == 'Load' else 'unloading')
    ax.axhline(0.0, color='0.6', lw=0.6)
    ax.set_xlim(0, DMAX)
    ax.set_xlabel('Displacement of top plate [mm]')
    ax.set_ylabel('Reaction force [N]')
    ax.set_title('%s - %s (up to %g mm)' % (prefix, label, DMAX))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(HERE, 'force_disp_%s.png' % prefix)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print('Wrote ' + out)
