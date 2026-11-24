"""Analyse, par type de lipide, des changements de feuillet (V4) et des desaccords avec la reference geometrique
(plan median local, grille ~3 nm). Usage :
.venv/bin/python scripts/analyse_types_trajectoire.py structure.gro traj.xtc topol.tpr sortie.txt [RES_EXCLUS]"""
import collections
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lfpy import *  # noqa
gro, traj, tpr, out = sys.argv[1:5]
excl = sys.argv[5].split(",") if len(sys.argv) > 5 else []
u = mda.Universe(gro, traj)
A = Analyzer(u, molnums=tpr_molnums(tpr), exclude=excl)
A.prepare()
h = A.head_atoms(auto=True)
T = u.atoms.resnames[h]
prev = None
chg = collections.Counter(); dchg = []
dis = collections.Counter(); dz_d = []
n = 0
for ts in u.trajectory:
    lab, _ = A.run(4)
    P = u.atoms.positions[h] / 10
    d = local_midplane_distance(P, ts.dimensions[:3] / 10)
    if prev is not None:
        if np.mean(lab == prev) < 0.5:
            lab = 1 - lab
        m = lab != prev
        for t in T[m]:
            chg[t] += 1
        dchg += list(np.abs(d[m]))
    prev = lab
    zr = (d > 0).astype(int)
    if np.mean(zr == lab) < 0.5:
        zr = 1 - zr
    bad = zr != lab
    for t in T[bad]:
        dis[t] += 1
    dz_d += list(np.abs(d[bad]))
    n += 1
L = ["Analyse par type de lipide : %s (%d trames, V4)" % (traj, n),
     "composition : %s" % dict(collections.Counter(T)),
     "changements de feuillet entre trames consecutives, moyenne par trame et par type : %s" %
     {k: round(v / max(1, n - 1), 2) for k, v in chg.items()},
     "  |distance au plan median local| des lipides qui changent (nm) : mediane %.2f, 90e centile %.2f" %
     ((np.median(dchg), np.percentile(dchg, 90)) if dchg else (0, 0)),
     "desaccords V4 / reference geometrique, moyenne par trame et par type : %s" %
     {k: round(v / n, 1) for k, v in dis.items()},
     "  |distance au plan median local| des desaccords (nm) : mediane %.2f, 90e centile %.2f, max %.2f" %
     (np.median(dz_d), np.percentile(dz_d, 90), np.max(dz_d))]
print("\n".join(L))
open(out, "w").write("\n".join(L) + "\n")
