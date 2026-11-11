# -*- coding: utf-8 -*-
"""
Multistable energy-trapping lattice (molded design of Shan et al., Adv. Mater. 2015)
Abaqus/CAE script -> 2D plane strain, Abaqus/Explicit quasi-static compression
followed by unloading, between frictional rigid plates (specimen not bonded).

Run:   abaqus cae noGUI=shan_multistable.py     (or CAE: File > Run Script)
Units: mm, N, tonne, s   (stress in MPa, density in t/mm^3)

Geometry is rebuilt parametrically from drawing Model-1a-Layout1 (106 x 108 mm,
4 beam layers x 8 tilted beams). Resulting beam: t = 1.01 mm, face length 6.29 mm,
inclination 66.1 deg from horizontal (t/L ~ 0.16).
"""

from abaqus import *
from abaqusConstants import *
from caeModules import *
import mesh

# =============================== USER PARAMETERS ===============================
MATERIAL  = 'ELASTIC50A'   # 'ELASTIC50A' or 'TPU90A'
DEPTH     = 30.0           # out-of-plane depth = plane-strain thickness [mm]
STROKE    = 30.0           # top-plate travel [mm]
T_LOAD    = 2.0            # loading step time [s]
T_UNLOAD  = 2.0            # unloading step time [s]
DT_TARGET = 2.0e-6         # mass-scaling target increment [s]; 0 -> no mass scaling
ALPHA     = 5.0            # mass-proportional Rayleigh damping [1/s]
K_OVER_MU = 100.0          # bulk/shear ratio (nu ~ 0.495). Shan et al. used 2500.
SIZE_BEAM = 0.2            # element length along the beams [mm]
N_THICK   = 5              # elements through the beam thickness (structured)
SIZE_GLOB = 0.8            # global seed [mm]
MU_PLATE  = 0.3            # friction coefficient, plates <-> specimen
MU_SELF   = 0.1            # friction coefficient, self-contact
NCPU      = 4
MODEL_NAME = JOB_NAME = 'Shan_Multistable'
SUBMIT    = False          # True -> submit and wait

# Hyperelastic constants (Yeoh: C10, C20, C30;  neo-Hookean: C10) [MPa]
MATERIALS = {
    # Formlabs Elastic 50A V2. 3-term Yeoh fitted to datasheet nominal stress
    # 0.9 / 1.7 / 3.4 MPa at 50 / 100 / 160 % strain. mu0 = 0.83 MPa.
    # Density not on datasheet: 1.1 g/cm^3 assumed.
    'ELASTIC50A': dict(rho=1.10e-9, C=(0.41373, 0.00782, 0.00339)),
    # TPU 90A (Formlabs SLS, XY). PROVISIONAL: neo-Hookean matched to 6.1 MPa at
    # 50 % strain only. TPU shows yield-like softening and ~20 % compression set,
    # which a hyperelastic model cannot capture. Replace with your coupon data.
    'TPU90A':     dict(rho=1.14e-9, C=(2.89,)),
}

# ================================== GEOMETRY ===================================
W, H      = 106.0, 108.0
LAYERS    = [(15.0, 24.5, 'up'), (34.5, 44.0, 'down'),     # (y_bottom, y_top, type)
             (64.0, 73.5, 'up'), (83.5, 93.0, 'down')]
HOLE_WIDE = 14.0     # wide end of each opening
HOLE_NARR = 8.9      # narrow end of each opening
SLANT_H   = 5.75     # vertical extent of the tilted faces
SLOT      = 5.0      # offset of the 'A' openings into the adjacent bar
A_CX      = (17.0, 41.0, 65.0, 89.0)   # openings sunk into the bar
B_CX      = (29.0, 53.0, 77.0)         # openings inside the beam layer
NOTCH_IN  = 12.0     # notch depth at the side faces (B opening centred at x = 5)


def opening(cx, y_narrow, y_wide):
    """Hexagonal opening: narrow end, tilted faces over SLANT_H, then straight sides."""
    sg = 1.0 if y_wide > y_narrow else -1.0
    ys = y_narrow + sg * SLANT_H
    hn, hw = HOLE_NARR / 2.0, HOLE_WIDE / 2.0
    return [(cx - hn, y_narrow), (cx + hn, y_narrow), (cx + hw, ys),
            (cx + hw, y_wide), (cx - hw, y_wide), (cx - hw, ys)]


def find_beams(segs, ext=0.05):
    """Pair parallel tilted faces into beams. For each beam return the two
    perpendicular cuts that bound its straight rectangular part (extended by `ext`
    into the openings), the rectangle centre, and mid-points of the face pieces."""
    beams = []
    for i in range(len(segs)):
        (ax, ay), (bx, by) = segs[i]
        L = ((bx - ax) ** 2 + (by - ay) ** 2) ** 0.5
        ux, uy = (bx - ax) / L, (by - ay) / L
        nx, ny = -uy, ux
        for j in range(i + 1, len(segs)):
            (cx, cy), (dx, dy) = segs[j]
            if abs((dx - cx) * uy - (dy - cy) * ux) > 1e-6:
                continue                                    # not parallel
            d = (cx - ax) * nx + (cy - ay) * ny
            if not 0.5 < abs(d) < 1.5:
                continue
            s1, s2 = sorted(((cx - ax) * ux + (cy - ay) * uy, (dx - ax) * ux + (dy - ay) * uy))
            lo, hi = max(0.0, s1), min(L, s2)
            if hi - lo < 1.0:
                continue
            pt = lambda s, off: (ax + s * ux + off * nx, ay + s * uy + off * ny)
            sg = 1.0 if d > 0 else -1.0
            cuts = [(pt(s, -sg * ext), pt(s, d + sg * ext)) for s in (lo, hi)]
            along, roots = [], []
            for off, ends in ((0.0, (0.0, L)), (d, (s1, s2))):
                along.append(pt(0.5 * (lo + hi), off))
                for e in ends:
                    if abs(e - lo) > 1e-6 and abs(e - hi) > 1e-6:
                        near = lo if abs(e - lo) < abs(e - hi) else hi
                        roots.append(pt(0.5 * (e + near), off))
            beams.append(dict(cuts=cuts, centre=pt(0.5 * (lo + hi), 0.5 * d),
                              cut_mids=[pt(s, 0.5 * d) for s in (lo, hi)],
                              along=along, roots=roots, length=hi - lo))
    return beams


def build_geometry():
    holes, slants = [], []
    for yb, yt, kind in LAYERS:
        for cx in A_CX:
            holes.append(opening(cx, yt - SLOT, yb - SLOT) if kind == 'up'
                         else opening(cx, yb + SLOT, yt + SLOT))
        for cx in B_CX:
            holes.append(opening(cx, yb, yt) if kind == 'up' else opening(cx, yt, yb))
    # right side face, bottom -> top (notch = B opening centred at x = W - 5)
    xn = W - (NOTCH_IN - (HOLE_WIDE - HOLE_NARR) / 2.0)   # 96.55
    xw = W - NOTCH_IN                                      # 94.0
    right = [(W, 0.0)]
    for yb, yt, kind in LAYERS:
        if kind == 'up':
            right += [(W, yb), (xn, yb), (xw, yb + SLANT_H), (xw, yt), (W, yt)]
        else:
            right += [(W, yb), (xw, yb), (xw, yt - SLANT_H), (xn, yt), (W, yt)]
    right.append((W, H))
    left = [(W - x, y) for x, y in reversed(right)]        # mirror, top -> bottom
    outer = right + left                                    # counter-clockwise
    for loop in [outer] + holes:
        n = len(loop)
        for i in range(n):
            (x1, y1), (x2, y2) = loop[i], loop[(i + 1) % n]
            if abs(x2 - x1) > 1e-9 and abs(y2 - y1) > 1e-9:
                slants.append(((x1, y1), (x2, y2)))
    return outer, holes, find_beams(slants)


# ================================ ABAQUS MODEL =================================

outer, holes, beams = build_geometry()
m = mdb.Model(name=MODEL_NAME, modelType=STANDARD_EXPLICIT)
if 'Model-1' in mdb.models.keys() and len(mdb.models.keys()) > 1:
    del mdb.models['Model-1']

# ---- lattice part
sk = m.ConstrainedSketch(name='lattice', sheetSize=400.0)
for loop in [outer] + holes:
    for i in range(len(loop)):
        sk.Line(point1=loop[i], point2=loop[(i + 1) % len(loop)])
p = m.Part(name='Lattice', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p.BaseShell(sketch=sk)

# partition each beam into a clean rectangle (structured mesh) + root transitions
tr = p.MakeSketchTransform(sketchPlane=p.faces[0], sketchPlaneSide=SIDE1,
                           origin=(0.0, 0.0, 0.0))
ps = m.ConstrainedSketch(name='beam_cuts', sheetSize=400.0, transform=tr)
for b in beams:
    for P, Q in b['cuts']:
        ps.Line(point1=P, point2=Q)
p.PartitionFaceBySketch(faces=p.faces, sketch=ps)

xyz = lambda pts: [((x, y, 0.0),) for x, y in pts]
p.Set(faces=p.faces, name='ALL')
beam_faces = p.faces.findAt(*xyz([b['centre'] for b in beams]))
along_edges = p.edges.findAt(*xyz([q for b in beams for q in b['along']]))
root_edges = p.edges.findAt(*xyz([q for b in beams for q in b['roots']]))
cut_edges = p.edges.findAt(*xyz([q for b in beams for q in b['cut_mids']]))
p.Set(faces=beam_faces, name='BEAMS')
p.Set(edges=along_edges, name='BEAM_FACES')

# ---- material and section
md = MATERIALS[MATERIAL]
mu0 = 2.0 * md['C'][0]
D1 = 2.0 / (K_OVER_MU * mu0)
mat = m.Material(name=MATERIAL)
mat.Density(table=((md['rho'],),))
if len(md['C']) == 3:
    mat.Hyperelastic(materialType=ISOTROPIC, testData=OFF, type=YEOH,
                     volumetricResponse=VOLUMETRIC_DATA,
                     table=(tuple(md['C']) + (D1, 0.0, 0.0),))
else:
    mat.Hyperelastic(materialType=ISOTROPIC, testData=OFF, type=NEO_HOOKE,
                     volumetricResponse=VOLUMETRIC_DATA, table=((md['C'][0], D1),))
if ALPHA > 0.0:
    mat.Damping(alpha=ALPHA)
m.HomogeneousSolidSection(name='SEC', material=MATERIAL, thickness=DEPTH)
p.SectionAssignment(region=p.sets['ALL'], sectionName='SEC')

# ---- mesh
n_along = max(1, int(round(beams[0]['length'] / SIZE_BEAM)))
p.seedPart(size=SIZE_GLOB, deviationFactor=0.1, minSizeFactor=0.1)
p.seedEdgeByNumber(edges=along_edges, number=n_along, constraint=FIXED)
p.seedEdgeByNumber(edges=cut_edges, number=N_THICK, constraint=FIXED)
p.seedEdgeBySize(edges=root_edges, size=SIZE_BEAM, deviationFactor=0.1, constraint=FINER)
p.setMeshControls(regions=p.faces, elemShape=QUAD, technique=FREE, algorithm=ADVANCING_FRONT)
p.setMeshControls(regions=beam_faces, elemShape=QUAD, technique=STRUCTURED)
p.setElementType(regions=(p.faces,), elemTypes=(
    mesh.ElemType(elemCode=CPE4R, elemLibrary=EXPLICIT, kinematicSplit=AVERAGE_STRAIN,
                  hourglassControl=ENHANCED, distortionControl=DEFAULT),
    mesh.ElemType(elemCode=CPE3, elemLibrary=EXPLICIT)))
p.generateMesh()
q = p.verifyMeshQuality(criterion=ANALYSIS_CHECKS)
n_bad = len(q['failedElements'])
print('Mesh check: %d elements, %d failed, %d warnings'
      % (len(p.elements), n_bad, len(q['warningElements'])))
if n_bad:
    p.Set(elements=q['failedElements'], name='BAD_ELEMENTS')
    raise RuntimeError('%d bad elements (part set BAD_ELEMENTS). Adjust SIZE_GLOB/SIZE_BEAM.'
                       % n_bad)


# ---- rigid plates (normal = left of the drawing direction -> both face the lattice)
def plate(name, y, x_from, x_to):
    s = m.ConstrainedSketch(name=name, sheetSize=400.0)
    s.Line(point1=(x_from, y), point2=(x_to, y))
    pl = m.Part(name=name, dimensionality=TWO_D_PLANAR, type=ANALYTIC_RIGID_SURFACE)
    pl.AnalyticRigidSurf2DPlanar(sketch=s)
    return pl, pl.ReferencePoint(point=(W / 2.0, y, 0.0)).id


pB, rpB = plate('PlateBot', 0.0, -20.0, W + 20.0)   # normal +y
pT, rpT = plate('PlateTop', H, W + 20.0, -20.0)     # normal -y

# ---- assembly
a = m.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
iL = a.Instance(name='Lattice-1', part=p, dependent=ON)
iB = a.Instance(name='PlateBot-1', part=pB, dependent=ON)
iT = a.Instance(name='PlateTop-1', part=pT, dependent=ON)
a.Set(referencePoints=(iB.referencePoints[rpB],), name='RP_BOT')
a.Set(referencePoints=(iT.referencePoints[rpT],), name='RP_TOP')
a.Surface(side1Edges=iB.edges, name='S_PLATE_BOT')
a.Surface(side1Edges=iT.edges, name='S_PLATE_TOP')
a.Surface(side1Edges=iL.edges.findAt(((W / 2.0, 0.0, 0.0),)), name='S_LAT_BOT')
a.Surface(side1Edges=iL.edges.findAt(((W / 2.0, H, 0.0),)), name='S_LAT_TOP')


# Self-contact: Explicit contact pairs in 2D need one continuous surface each, so
# each closed boundary loop (outer + 28 openings) gets its own surface. Contact can
# only happen across a void, i.e. between edges of the same loop, so this is complete.
def on_loop(pt, loop, tol=1e-4):
    for i in range(len(loop)):
        (x1, y1), (x2, y2) = loop[i], loop[(i + 1) % len(loop)]
        dx, dy = x2 - x1, y2 - y1
        s = max(0.0, min(1.0, ((pt[0] - x1) * dx + (pt[1] - y1) * dy) / (dx * dx + dy * dy)))
        if ((pt[0] - x1 - s * dx) ** 2 + (pt[1] - y1 - s * dy) ** 2) ** 0.5 < tol:
            return True
    return False


loops = [outer] + holes
groups = [[] for _ in loops]
for e in iL.edges:
    if len(e.getFaces()) != 1:
        continue                          # partition line (internal edge)
    pt = e.pointOn[0]
    k = [i for i, lp in enumerate(loops) if on_loop(pt, lp)]
    if len(k) != 1:
        raise RuntimeError('Boundary edge at %s not matched to a loop' % (pt,))
    groups[k[0]].append(pt)
self_surfaces = []
for k, pts in enumerate(groups):
    name = 'S_LOOP_%02d' % k
    a.Surface(side1Edges=iL.edges.findAt(*[(q,) for q in pts]), name=name)
    self_surfaces.append(name)

# ---- steps
ms = PREVIOUS_STEP
if DT_TARGET > 0.0:
    ms = ((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 0.0, DT_TARGET, BELOW_MIN,
           0, 0, 0.0, 0.0, 0, None),)
m.ExplicitDynamicsStep(name='Load', previous='Initial', timePeriod=T_LOAD, massScaling=ms)
m.ExplicitDynamicsStep(name='Unload', previous='Load', timePeriod=T_UNLOAD)
m.SmoothStepAmplitude(name='AMP_LOAD', timeSpan=STEP, data=((0.0, 0.0), (T_LOAD, 1.0)))
m.SmoothStepAmplitude(name='AMP_UNLOAD', timeSpan=STEP, data=((0.0, 1.0), (T_UNLOAD, 0.0)))

# ---- boundary conditions
m.EncastreBC(name='FIX_BOT', createStepName='Initial', region=a.sets['RP_BOT'])
m.DisplacementBC(name='DRIVE_TOP', createStepName='Initial', region=a.sets['RP_TOP'],
                 u1=SET, u2=SET, ur3=SET)
bc = m.boundaryConditions['DRIVE_TOP']
bc.setValuesInStep(stepName='Load', u2=-STROKE, amplitude='AMP_LOAD')
bc.setValuesInStep(stepName='Unload', u2=-STROKE, amplitude='AMP_UNLOAD')

# ---- contact (plates not bonded: lattice can lift off, as in the unattached tests)
for name, mu in (('IP_PLATE', MU_PLATE), ('IP_SELF', MU_SELF)):
    ip = m.ContactProperty(name)
    ip.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON,
                      constraintEnforcementMethod=DEFAULT)
    ip.TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC,
                          slipRateDependency=OFF, pressureDependency=OFF,
                          temperatureDependency=OFF, dependencies=0, table=((mu,),),
                          shearStressLimit=None, maximumElasticSlip=FRACTION,
                          fraction=0.005, elasticSlipStiffness=None)


def plate_contact(name, main, secondary):
    kw = dict(name=name, createStepName='Initial', sliding=FINITE,
              mechanicalConstraint=KINEMATIC, interactionProperty='IP_PLATE')
    try:   # Abaqus 2022+
        m.SurfaceToSurfaceContactExp(main=main, secondary=secondary, **kw)
    except TypeError:   # older releases
        m.SurfaceToSurfaceContactExp(master=main, slave=secondary, **kw)


plate_contact('C_BOT', a.surfaces['S_PLATE_BOT'], a.surfaces['S_LAT_BOT'])
plate_contact('C_TOP', a.surfaces['S_PLATE_TOP'], a.surfaces['S_LAT_TOP'])
for name in self_surfaces:
    m.SelfContactExp(name='C_SELF_' + name[-2:], createStepName='Initial',
                     surface=a.surfaces[name], mechanicalConstraint=PENALTY,
                     interactionProperty='IP_SELF')

# ---- output
m.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'LE', 'U', 'V', 'RF', 'CSTRESS'), numIntervals=100)
m.historyOutputRequests['H-Output-1'].setValues(
    variables=('ALLKE', 'ALLIE', 'ALLSE', 'ALLVD', 'ALLAE', 'ALLFD', 'ALLWK', 'ETOTAL'),
    numIntervals=1000)
m.HistoryOutputRequest(name='H-RP_TOP', createStepName='Load', variables=('U2', 'RF2'),
                       region=a.sets['RP_TOP'], numIntervals=1000)

# ---- job
mdb.Job(name=JOB_NAME, model=MODEL_NAME, type=ANALYSIS, explicitPrecision=DOUBLE_PLUS_PACK,
        nodalOutputPrecision=FULL, numCpus=NCPU, numDomains=NCPU,
        multiprocessingMode=DEFAULT)
mdb.saveAs(pathName=JOB_NAME + '.cae')
mdb.jobs[JOB_NAME].writeInput(consistencyChecking=OFF)
print('Elements: %d   material: %s   mu0 = %.3f MPa   D1 = %.4g'
      % (len(p.elements), MATERIAL, mu0, D1))
if SUBMIT:
    mdb.jobs[JOB_NAME].submit(consistencyChecking=OFF)
    mdb.jobs[JOB_NAME].waitForCompletion()
