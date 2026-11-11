# -*- coding: utf-8 -*-
"""
Force-displacement post-processing for the tilted-beam lattice mesh study.

Reads U2, RF2 at RP_TOP (history output) for steps 'Load' and 'Unload' from every
ShanMesh*.odb and plots loading (solid) and unloading (dashed) of all meshes
on one graph. Displacement is compression (down +), force is compressive +.

Run:   python force_disp.py [odb ...]        (default: every *.odb in this folder)
  - extraction needs odbAccess, so the script calls `abaqus python` on itself
    for any ODB whose <job>_fd.csv is missing or older than the ODB;
  - plotting is done with the normal Python (matplotlib).
Writes:  <job>_fd.csv   (step, disp_mm, force_N)   and   force_disp.png
"""
import glob
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = ('Load', 'Unload')
RP_KEY = 'PLATETOP'                  # history region "Node PLATETOP-1.1" (RP_TOP)
ABAQUS = shutil.which('abaqus') or r'C:\Program Files\SIMULIA\Commands\abaqus.bat'
NAME_RE = re.compile(r'ShanMesh([\d._]+)mm([\d._]+)mm$')


def csv_path(odb):
    return os.path.splitext(odb)[0] + '_fd.csv'


# ---------------------------------------------------------------- extraction (abaqus python)
def extract(odb_path):
    from odbAccess import openOdb
    odb = openOdb(odb_path, readOnly=True)
    with open(csv_path(odb_path), 'w') as fh:
        fh.write('step,disp_mm,force_N\n')
        for step in STEPS:
            if step not in odb.steps.keys():
                continue
            regions = odb.steps[step].historyRegions
            keys = [k for k in regions.keys()
                    if {'U2', 'RF2'} <= set(regions[k].historyOutputs.keys())]
            if not keys:
                raise RuntimeError('No U2/RF2 history output in step %s of %s' % (step, odb_path))
            ho = regions[([k for k in keys if RP_KEY in k.upper()] or keys)[0]].historyOutputs
            for (_, u), (_, f) in zip(ho['U2'].data, ho['RF2'].data):
                fh.write('%s,%.6g,%.6g\n' % (step, -u, -f))
    odb.close()
    print('Wrote ' + csv_path(odb_path))


# ---------------------------------------------------------------- plotting (normal python)
def mesh_size(job):
    """'ShanMesh5mm0_5mm' -> (5.0, 0.5): (global seed, beam seed) in mm, or None."""
    m = NAME_RE.match(job)
    return tuple(float(g.replace('_', '.')) for g in m.groups()) if m else None


def read_csv(fname):
    out = {}
    with open(fname) as fh:
        next(fh)
        for line in fh:
            s, u, f = line.strip().split(',')
            d = out.setdefault(s, ([], []))
            d[0].append(float(u))
            d[1].append(float(f))
    return out


def plot(odbs, out):
    import numpy as np
    import matplotlib.pyplot as plt

    jobs = [os.path.splitext(os.path.basename(o))[0] for o in odbs]
    # coarse -> fine; names that don't parse go last, alphabetically
    order = sorted(range(len(odbs)), key=lambda i: (
        (0, -mesh_size(jobs[i])[1], -mesh_size(jobs[i])[0]) if mesh_size(jobs[i])
        else (1, jobs[i])))
    colors = plt.cm.viridis(np.linspace(0.85, 0.0, len(odbs)))     # finer = darker

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for c, i in zip(colors, order):
        data = read_csv(csv_path(odbs[i]))
        sz = mesh_size(jobs[i])
        lab = 'global %g mm, beam %g mm' % sz if sz else jobs[i]
        if 'Load' in data:
            ax.plot(*data['Load'], color=c, lw=1.2, label=lab)
        if 'Unload' in data:
            ax.plot(*data['Unload'], color=c, lw=1.2, ls='--')
    ax.plot([], [], color='0.3', lw=1.2, label='loading')
    ax.plot([], [], color='0.3', lw=1.2, ls='--', label='unloading')
    ax.axhline(0.0, color='0.6', lw=0.6)
    ax.set_xlabel('Displacement of top plate [mm]')
    ax.set_ylabel('Reaction force [N]')
    ax.set_title('Force-displacement, mesh convergence')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print('Wrote ' + out)


def main(args):
    odbs = [os.path.abspath(a) for a in args] or sorted(glob.glob(os.path.join(HERE, '*.odb')))
    if not odbs:
        sys.exit('No .odb files found.')
    stale = [o for o in odbs if not os.path.exists(csv_path(o))
             or os.path.getmtime(csv_path(o)) < os.path.getmtime(o)]
    if stale:
        subprocess.check_call([ABAQUS, 'python', os.path.abspath(__file__), '--extract'] + stale,
                              cwd=HERE)
    plot(odbs, os.path.join(HERE, 'force_disp.png'))


if __name__ == '__main__':
    if '--extract' in sys.argv:                  # inside `abaqus python`
        for p in sys.argv[sys.argv.index('--extract') + 1:]:
            extract(p)
    else:
        main(sys.argv[1:])
