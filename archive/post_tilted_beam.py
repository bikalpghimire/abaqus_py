# -*- coding: utf-8 -*-
"""
Extract results of the tilted-beam lattice model (ShanMesh*.inp) from its ODB.

Model (from the .inp): 2D lattice (CPE4R, 30 mm deep) between two rigid plates.
RP_TOP is driven down 38 mm in step 'Load' and back to 0 in step 'Unload'
(Abaqus/Explicit). History output: RF2, U2 at RP_TOP; whole-model energies.

Run:   abaqus python post_tilted_beam.py [<job>.odb]     (default: ODB_PATH)
Writes next to the ODB:
    <job>_load.csv       time, disp_mm, force_N      (loading)
    <job>_unload.csv     time, disp_mm, force_N      (unloading, if present)
    <job>_energy.csv     ALLSE, ALLKE, ALLIE, ALLVD, ALLFD, ... vs time, all steps
    <job>_summary.json   key results, read by plot_tilted_beam.py
Displacement is compression of the lattice (downward +), force is compressive +.
Units follow the model: mm, N, tonne, s  ->  energy in N*mm = mJ.
"""
import json
import os
import sys
from odbAccess import openOdb

ODB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'Shan_Mesh5mm_0_5mm.odb')
RP_KEY = 'PLATETOP'                  # history region of RP_TOP ("Node PLATETOP-1.1")
LOAD_STEP, UNLOAD_STEP = 'Load', 'Unload'
ENERGY_KEYS = ['ALLWK', 'ALLIE', 'ALLSE', 'ALLKE', 'ALLVD', 'ALLFD', 'ALLAE',
               'ALLPD', 'ALLSD', 'ETOTAL']


def path_work(u, f):
    """Work of f along the path u (trapezoidal); sign follows the direction of travel."""
    return sum(0.5 * (f[i] + f[i - 1]) * (u[i] - u[i - 1]) for i in range(1, len(u)))


def zero_crossings(u, f):
    """Displacements where f changes sign, with the sign it changes to."""
    out = []
    for i in range(1, len(f)):
        if f[i - 1] > 0.0 >= f[i] or f[i - 1] < 0.0 <= f[i]:
            s = f[i - 1] / (f[i - 1] - f[i])
            out.append((u[i - 1] + s * (u[i] - u[i - 1]), -1 if f[i] <= 0.0 else 1))
    return out


def rp_curve(odb, step):
    """Total time, displacement (down +) and force (compression +) at RP_TOP."""
    regions = odb.steps[step].historyRegions
    cands = [k for k in regions.keys()
             if 'U2' in regions[k].historyOutputs.keys()
             and 'RF2' in regions[k].historyOutputs.keys()]
    if not cands:
        raise RuntimeError('No U2/RF2 history output in step %s' % step)
    key = ([k for k in cands if RP_KEY in k.upper()] or cands)[0]
    ho = regions[key].historyOutputs
    du, df = ho['U2'].data, ho['RF2'].data
    n = min(len(du), len(df))
    t0 = odb.steps[step].totalTime
    return ([t0 + du[i][0] for i in range(n)],
            [-du[i][1] for i in range(n)],
            [-df[i][1] for i in range(n)])


def write_curve(fname, t, u, f):
    with open(fname, 'w') as fh:
        fh.write('time,disp_mm,force_N\n')
        for row in zip(t, u, f):
            fh.write('%.6g,%.6g,%.6g\n' % row)


def energy_histories(odb, steps):
    """{name: [(total_time, value), ...]} over all steps (Explicit energies are cumulative)."""
    out = {}
    for s in steps:
        hr = odb.steps[s].historyRegions
        if 'Assembly ASSEMBLY' not in hr.keys():
            continue
        ho = hr['Assembly ASSEMBLY'].historyOutputs
        t0 = odb.steps[s].totalTime
        for k in ENERGY_KEYS:
            if k in ho.keys():
                out.setdefault(k, []).extend((t0 + t, v) for t, v in ho[k].data)
    return out


def write_energy(fname, odb, steps):
    rows, cols = [], []
    for s in steps:
        hr = odb.steps[s].historyRegions
        if 'Assembly ASSEMBLY' not in hr.keys():
            continue
        ho = hr['Assembly ASSEMBLY'].historyOutputs
        data = dict((k, ho[k].data) for k in ENERGY_KEYS if k in ho.keys())
        if not data:
            continue
        cols += [k for k in ENERGY_KEYS if k in data and k not in cols]
        m = min(len(d) for d in data.values())
        tt = list(data.values())[0]
        t0 = odb.steps[s].totalTime
        for i in range(m):
            rows.append((s, t0 + tt[i][0], dict((k, data[k][i][1]) for k in data)))
    if not rows:
        return
    with open(fname, 'w') as fh:
        fh.write('step,time,' + ','.join(cols) + '\n')
        for s, t, d in rows:
            fh.write('%s,%.6g,%s\n' % (s, t, ','.join(
                '%.6g' % d[k] if k in d else '' for k in cols)))


def analyze(path):
    folder = os.path.dirname(os.path.abspath(path))
    job = os.path.splitext(os.path.basename(path))[0]
    base = os.path.join(folder, job)

    odb = openOdb(path, readOnly=True)
    steps = list(odb.steps.keys())
    explicit = 'EXPLICIT' in str(odb.jobData.analysisCode).upper()
    load = LOAD_STEP if LOAD_STEP in steps else steps[0]
    unload = UNLOAD_STEP if UNLOAD_STEP in steps else None

    # ---- loading ----
    t, u, f = rp_curve(odb, load)
    write_curve(base + '_load.csv', t, u, f)
    i_pk = max(range(len(f)), key=lambda i: f[i])
    e_load = path_work(u, f)

    # snap-through / negative force during loading (bistable-type behaviour)
    zc = [z for z in zero_crossings(u, f) if z[0] > 0.0]
    neg = [z[0] for z in zc if z[1] < 0]
    snap = min(f[i_pk:]) < 0.95 * f[i_pk] and i_pk < len(f) - 1

    # ---- unloading ----
    S = dict(job=job, solver='EXPLICIT' if explicit else 'STANDARD',
             load_step=load, unload_step=unload,
             u_peak=u[i_pk], f_peak=f[i_pk], u_max=max(u), f_at_umax=f[u.index(max(u))],
             e_load=e_load, e_return=None, e_diss=None, eta=None,
             f_end_unload=None, u_zero_force=neg[0] if neg else None,
             snap_through=bool(snap))
    if unload:
        tb, ub, fb = rp_curve(odb, unload)
        write_curve(base + '_unload.csv', tb, ub, fb)
        e_return = -path_work(ub, fb)               # work given back by the lattice
        S.update(e_return=e_return, e_diss=e_load - e_return,
                 eta=(e_load - e_return) / e_load if e_load else None,
                 f_end_unload=fb[-1])

    # ---- energies ----
    write_energy(base + '_energy.csv', odb, steps)
    E = energy_histories(odb, steps)
    last = dict((k, v[-1][1]) for k, v in E.items())
    S['energy_end'] = last

    # ---- report ----
    print('')
    print('%s   (%s; steps: %s)' % (job, S['solver'].title(), ', '.join(steps)))
    print('Max compression     u_max  = %.4g mm   (P there = %.4g N)' % (S['u_max'], S['f_at_umax']))
    print('Peak force          P_peak = %.4g N  at u = %.4g mm' % (f[i_pk], u[i_pk]))
    if snap:
        print('  force drops >5 % after the peak: snap-through / instability')
    if neg:
        print('  force turns negative at u = %.4g mm (lattice pulls the plate: bistable cells)' % neg[0])
    print('Energy absorbed on loading   E_load   = %.4g mJ' % e_load)
    if unload:
        print('Energy returned on unloading E_return = %.4g mJ' % S['e_return'])
        print('Dissipated (loop area)       E_diss   = %.4g mJ   (eta = E_diss/E_load = %.1f %%)'
              % (S['e_diss'], 100.0 * S['eta']))
        print('Force at end of unloading = %.4g N' % S['f_end_unload'])
        if S['e_diss'] > 0:
            for k, lab in [('ALLVD', '*Damping alpha + bulk viscosity'),
                           ('ALLFD', 'friction, self-contact')]:
                if k in last:
                    print('  %s = %.4g mJ = %5.1f %% of E_diss   (%s)'
                          % (k, last[k], 100.0 * last[k] / S['e_diss'], lab))

    # ---- numerical-quality checks ----
    if 'ALLIE' in E:
        ie_max = max(v for _, v in E['ALLIE'])
        pairs = lambda a: [(x[1], y[1]) for x, y in zip(E[a], E['ALLIE']) if y[1] > 0.01 * ie_max]
        if explicit and 'ALLKE' in E:
            r = max(a / b for a, b in pairs('ALLKE'))
            print('Max ALLKE / ALLIE = %.1f %%   (quasi-static if < ~5 %%; snaps spike briefly)' % (100 * r))
        if 'ALLAE' in E:
            r = max(a / b for a, b in pairs('ALLAE'))
            print('Max ALLAE / ALLIE = %.1f %%   (hourglass control; keep < ~5 %%)' % (100 * r))
    if 'ETOTAL' in last and 'ALLWK' in last and last['ALLWK']:
        print('ETOTAL / ALLWK at end = %.2f %%   (energy balance; should stay ~0)'
              % (100 * last['ETOTAL'] / max(abs(v) for _, v in E['ALLWK'])))

    with open(base + '_summary.json', 'w') as fh:
        json.dump(S, fh, indent=2)
    print('Wrote %s_load.csv%s, _energy.csv, _summary.json'
          % (job, ', _unload.csv' if unload else ''))
    odb.close()
    return S


if __name__ == '__main__':
    analyze(sys.argv[1] if len(sys.argv) > 1 else ODB_PATH)
