"""
run_trajectory.py -- point 5 et 6 de l'enonce : application de la version finale (V4)
a une trajectoire complete, mesure du temps de calcul, comparaison avec V1,
LeafletFinder (MDAnalysis) et FATSLiM.

Usage :
  .venv/bin/python scripts/run_trajectory.py [--gro F.gro] [--traj F.xtc] [--repeat R]
Par defaut : trajectoire/dppc_vesicle.gro + trajectoire/dppc_vesicle.xtc
Pour la trajectoire officielle du cours : --gro <structure.gro> --traj <trajectoire.xtc>

Sorties (resultats/trajectoire/) :
  par_trame.csv           temps et accords par trame et par methode
  feuillets_V4.txt        feuillets de chaque trame (V4), format leaflets.txt
  resume.txt              synthese (temps total, debit, extrapolation a 100 000 trames)
"""
import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lfpy import RESULTS, ROOT, Analyzer, accuracy  # noqa: E402
import MDAnalysis as mda  # noqa: E402
from MDAnalysis.analysis.leaflet import LeafletFinder  # noqa: E402
from benchmark_structures import run_fatslim, labels_from_groups  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--gro", default=os.path.join(ROOT, "trajectoire", "dppc_vesicle.gro"))
ap.add_argument("--traj", default=os.path.join(ROOT, "trajectoire", "dppc_vesicle.xtc"))
ap.add_argument("--repeat", type=int, default=200, help="passes sur les trames en memoire pour mesurer le debit")
ap.add_argument("--ref", default=None, help="fichier leaflets.txt de reference pour la trame 0 (optionnel)")
ap.add_argument("--out", default="trajectoire", help="sous-dossier de resultats/")
ap.add_argument("--tpr", default=None, help="fichier TPR : decoupage exact en molecules (lipides sur plusieurs residus)")
ap.add_argument("--exclude", default="", help="residus a exclure, separes par des virgules (p. ex. S14)")
ap.add_argument("--zref", action="store_true",
                help="bicouche plane : reference geometrique (tete au-dessus / au-dessous du plan median)")
ap.add_argument("--frames", default=None, help="sous-ensemble de trames debut:fin:pas (p. ex. 0:1001:1)")
args = ap.parse_args()
if args.ref is None and os.path.basename(args.gro) == "dppc_vesicle.gro":
    from lfpy import DATA
    args.ref = os.path.join(DATA, "dppc_vesicle", "leaflets.txt")
OUT = os.path.join(RESULTS, args.out)
os.makedirs(OUT, exist_ok=True)
NCPU = os.cpu_count()

u = mda.Universe(args.gro, args.traj)
MOL = None
if args.tpr:
    from lfpy import tpr_molnums
    MOL = tpr_molnums(args.tpr)
EXCL = [x for x in args.exclude.split(",") if x]
mk = lambda: Analyzer(u, nthreads=1, molnums=MOL, exclude=EXCL)
SL = slice(*[int(x) if x else None for x in args.frames.split(":")]) if args.frames else slice(None)
A = mk()
A1 = mk()
nf = len(u.trajectory[SL])
print("Trajectoire : %s (%d trames, %d atomes, %d lipides)" % (args.traj, nf, len(u.atoms), A.nlip))

# ---- 1) lecture des trames (I/O, mesure separee)
t0 = time.perf_counter()
frames = []
times_ps = []
for ts in u.trajectory[SL]:
    frames.append(Analyzer.frame_arrays(ts))
    times_ps.append(ts.time)
t_read = time.perf_counter() - t0

# ---- 2) detection des atomes de reference (une seule fois, sur la 1re trame)
t_prep = A.prepare(*frames[0])
A1.prepare(*frames[0])

# ---- 3) V4 et V1 trame par trame
labs4, t4, labs1, t1, cut1 = [], [], [], [], []
for xyz, box in frames:
    lab, info = A.run(4, xyz=xyz, box=box)
    labs4.append(lab); t4.append(info["time"])
    lab, info = A1.run(1, xyz=xyz, box=box)
    labs1.append(lab); t1.append(info["time"]); cut1.append(info["cutoff"])

# coherence des numeros de feuillets d'une trame a l'autre (feuillet 0 = le plus grand a la trame 0)
def align(lab, ref):
    if np.mean(lab[(lab >= 0) & (ref >= 0)] == ref[(lab >= 0) & (ref >= 0)]) < 0.5:
        return np.where(lab >= 0, 1 - lab, lab)
    return lab

labs4 = [align(l, labs4[0]) for l in labs4]
labs1 = [align(l, labs4[0]) for l in labs1]

# ---- 4) LeafletFinder (seuil 15 A, PBC) trame par trame
heads = A.head_atoms(auto=False)
head2lip = np.full(len(u.atoms), -1); head2lip[heads] = np.arange(A.nlip)
ag = u.atoms[heads]
labsLF, tLF = [], []
for k, ts in enumerate(u.trajectory[SL]):
    t0 = time.perf_counter()
    lf = LeafletFinder(u, ag, cutoff=15.0, pbc=True)
    groups = [head2lip[g.ix] for g in lf.groups()]
    tLF.append(time.perf_counter() - t0)
    labsLF.append(align(labels_from_groups(groups, A.nlip), labs4[0]))

# ---- 5) FATSLiM sur toute la trajectoire (1 thread et tous les coeurs)
hd = {int(h): i for i, h in enumerate(heads)}
fat = {}
for nt in (1, NCPU):
    labs, fms, wall, log = run_fatslim(args.gro, heads, hd, nt, A.nlip, traj=args.traj)
    sel = list(range(len(u.trajectory)))[SL]
    fat[nt] = ([labs[k] for k in sel if k < len(labs)], [fms[k] for k in sel if k < len(fms)], wall)
labsF = [align(l, labs4[0]) for l in fat[1][0]]

# ---- 6) debit soutenu : R passes sur les trames en memoire, sequentiel puis trames en parallele
R = max(1, min(args.repeat, 2200 // nf))  # environ 2000 analyses au total
seq = frames * R
t0 = time.perf_counter()
for xyz, box in seq:
    A.run(4, xyz=xyz, box=box)
t_seq = time.perf_counter() - t0

pool_an = [mk() for _ in range(NCPU)]
for P in pool_an:
    P.prepare(*frames[0])
chunks = [seq[i::NCPU] for i in range(NCPU)]

def work(i):
    for xyz, box in chunks[i]:
        pool_an[i].run(4, xyz=xyz, box=box)

t0 = time.perf_counter()
with ThreadPoolExecutor(NCPU) as ex:
    list(ex.map(work, range(NCPU)))
t_par = time.perf_counter() - t0

# ---- sorties
ref = None
if args.ref and os.path.exists(args.ref):
    from lfpy import read_groups
    pos = {m: i for i, m in enumerate(A.molids)}
    ref = np.full(A.nlip, -1)
    for k, (_, ids) in enumerate(read_groups(args.ref)):
        for x in ids:
            if x in pos:
                ref[pos[x]] = k

# reference geometrique pour une bicouche (quasi) plane : la tete (atome nomme) est-elle au-dessus ou au-dessous
# du plan median LOCAL ? Le plan median est estime sur une grille d'environ 3 nm x 3 nm (suit les ondulations).
zrefs = None
if args.zref:
    from lfpy import local_midplane_distance
    zrefs = []
    for k, (xyz, box) in enumerate(frames):
        d = local_midplane_distance(xyz[heads], np.diag(box))
        zrefs.append((d > 0).astype(int))
    # meme numerotation que V4 a la trame 0
    if np.mean(zrefs[0] == labs4[0]) < 0.5:
        zrefs = [1 - r for r in zrefs]

def agree(a, b):
    m = (a >= 0) | (b >= 0)
    return np.mean(a[m] == b[m])

with open(os.path.join(OUT, "par_trame.csv"), "w") as f:
    f.write("trame,temps_ps,V4_feuillet1,V4_feuillet2,V4_non_assignes,V4_ms,V1_ms,V1_seuil_nm,LeafletFinder_ms,"
            "FATSLiM_ms_1thr,FATSLiM_ms_%dthr,accord_V4_FATSLiM,accord_V4_LeafletFinder,accord_V4_V1,"
            "changements_V4_vs_trame_prec,accord_V4_reference_trame0,zref_V4,zref_V1,zref_LeafletFinder,zref_FATSLiM\n" % NCPU)
    for k in range(nf):
        l = labs4[k]
        chg = int(np.sum(l != labs4[k - 1])) if k else 0
        accr = accuracy(l, ref) if (ref is not None and k == 0) else ""
        zr = ["", "", "", ""]
        if zrefs is not None:
            zr = ["%.6f" % accuracy(L_[k], zrefs[k]) if k < len(L_) else "" for L_ in (labs4, labs1, labsLF, labsF)]
        f.write("%d,%.1f,%d,%d,%d,%.4f,%.4f,%.3f,%.3f,%s,%s,%.6f,%.6f,%.6f,%d,%s,%s\n" % (
            k, times_ps[k], np.sum(l == 0), np.sum(l == 1), np.sum(l < 0), t4[k] * 1e3, t1[k] * 1e3, cut1[k],
            tLF[k] * 1e3, fat[1][1][k] if k < len(fat[1][1]) else "", fat[NCPU][1][k] if k < len(fat[NCPU][1]) else "",
            agree(l, labsF[k]) if k < len(labsF) else np.nan, agree(l, labsLF[k]), agree(l, labs1[k]), chg,
            ("%.6f" % accr) if accr != "" else "", ",".join(zr)))

with open(os.path.join(OUT, "feuillets_V4.txt"), "w") as f:
    for k in range(nf):
        for L in sorted(set(labs4[k].tolist())):
            ids = A.molids[labs4[k] == L]
            name = "non_assignes" if L < 0 else "leaflet_%d" % (L + 1)
            f.write("[ trame_%04d_%s ]\n" % (k, name))
            for i in range(0, len(ids), 15):
                f.write("".join("%6d" % x for x in ids[i:i + 15]) + "\n")

n_tot = len(seq)
lines = [
    "Trajectoire               : %s" % args.traj,
    "Trames / atomes / lipides : %d / %d / %d" % (nf, len(u.atoms), A.nlip),
    "Lecture des trames (I/O, MDAnalysis) : %.1f ms (%.2f ms/trame)" % (t_read * 1e3, t_read * 1e3 / nf),
    "Detection des atomes de reference (1 fois) : %.1f ms" % (t_prep * 1e3),
    "",
    "Temps d'analyse sur l'ensemble de la trajectoire (%d trames) :" % nf,
    "  V4 (version finale, 1 thread)     : total %.2f ms, %.3f ms/trame" % (sum(t4) * 1e3, np.mean(t4) * 1e3),
    "  V1 (graphe + seuil auto, 1 thread): total %.2f ms, %.3f ms/trame" % (sum(t1) * 1e3, np.mean(t1) * 1e3),
    "  LeafletFinder (15 A, PBC)         : total %.1f ms, %.2f ms/trame" % (sum(tLF) * 1e3, np.mean(tLF) * 1e3),
    "  FATSLiM 1 thread                  : total %.0f ms (somme par trame), processus %.0f ms" % (sum(fat[1][1]), fat[1][2] * 1e3),
    "  FATSLiM %d threads                 : total %.0f ms (somme par trame), processus %.0f ms" % (NCPU, sum(fat[NCPU][1]), fat[NCPU][2] * 1e3),
    "",
    "Debit soutenu V4 (%d trames = %d passes, trames en memoire) :" % (n_tot, R),
    "  sequentiel 1 thread      : %.3f ms/trame -> %.0f trames/s" % (t_seq * 1e3 / n_tot, n_tot / t_seq),
    "  %d trames en parallele   : %.3f ms/trame -> %.0f trames/s" % (NCPU, t_par * 1e3 / n_tot, n_tot / t_par),
    "  extrapolation 100 000 trames : %.1f s (1 thread), %.1f s (%d threads), hors lecture du fichier" % (
        t_seq / n_tot * 1e5, t_par / n_tot * 1e5, NCPU),
    "",
    "Fiabilite :",
    "  accord V4 / reference (trame 0)       : %s" % (("%.4f" % accuracy(labs4[0], ref)) if ref is not None else "n/a"),
    "  accord moyen V4 / FATSLiM             : %.4f (min %.4f)" % (np.mean([agree(a, b) for a, b in zip(labs4, labsF)]),
                                                                 np.min([agree(a, b) for a, b in zip(labs4, labsF)])),
    "  accord moyen V4 / LeafletFinder       : %.4f" % np.mean([agree(a, b) for a, b in zip(labs4, labsLF)]),
    "  accord moyen V4 / V1                  : %.4f" % np.mean([agree(a, b) for a, b in zip(labs4, labs1)]),
    "  lipides changeant de feuillet entre trames (V4) : %s" % [int(np.sum(labs4[k] != labs4[k - 1])) for k in range(1, nf)],
    "  tailles des feuillets (V4) : %s" % [(int(np.sum(l == 0)), int(np.sum(l == 1)), int(np.sum(l < 0))) for l in labs4],
]
ZACC = {}
if zrefs is not None:
    for name, L_ in (("V4", labs4), ("V1", labs1), ("LeafletFinder", labsLF), ("FATSLiM", labsF)):
        a = [accuracy(L_[k], zrefs[k]) for k in range(min(nf, len(L_)))]
        ZACC[name] = (float(np.mean(a)), float(np.min(a)))
    lines.append("  accord avec la reference geometrique (tete au-dessus/au-dessous du plan median), moyenne (min) :")
    for name, (m_, mn_) in ZACC.items():
        lines.append("    %-14s : %.4f (%.4f)" % (name, m_, mn_))
txt = "\n".join(lines)
print(txt)
open(os.path.join(OUT, "resume.txt"), "w").write(txt + "\n")

import json
summary = dict(
    trajectoire=os.path.relpath(args.traj, ROOT), structure=os.path.relpath(args.gro, ROOT), n_trames=nf,
    n_atomes=len(u.atoms), n_lipides=int(A.nlip), duree_ps=float(times_ps[-1] - times_ps[0]),
    lecture_ms=t_read * 1e3, detection_tetes_ms=t_prep * 1e3,
    V4_total_ms=sum(t4) * 1e3, V4_ms_trame=float(np.mean(t4) * 1e3),
    V1_total_ms=sum(t1) * 1e3, V1_ms_trame=float(np.mean(t1) * 1e3),
    LF_total_ms=sum(tLF) * 1e3, LF_ms_trame=float(np.mean(tLF) * 1e3),
    FATSLiM1_total_ms=float(sum(fat[1][1])), FATSLiM1_processus_ms=fat[1][2] * 1e3,
    FATSLiMN_total_ms=float(sum(fat[NCPU][1])), FATSLiMN_processus_ms=fat[NCPU][2] * 1e3, ncpu=NCPU,
    debit_trames=n_tot, debit_seq_ms_trame=t_seq * 1e3 / n_tot, debit_par_ms_trame=t_par * 1e3 / n_tot,
    extrapol_100k_seq_s=t_seq / n_tot * 1e5, extrapol_100k_par_s=t_par / n_tot * 1e5,
    accord_ref_trame0=(float(accuracy(labs4[0], ref)) if ref is not None else None),
    accord_V4_FATSLiM_moy=float(np.mean([agree(a, b) for a, b in zip(labs4, labsF)])),
    accord_V4_FATSLiM_min=float(np.min([agree(a, b) for a, b in zip(labs4, labsF)])),
    accord_V4_LF_moy=float(np.mean([agree(a, b) for a, b in zip(labs4, labsLF)])),
    accord_V4_V1_moy=float(np.mean([agree(a, b) for a, b in zip(labs4, labs1)])),
    changements=[int(np.sum(labs4[k] != labs4[k - 1])) for k in range(1, nf)],
    tailles=[(int(np.sum(l == 0)), int(np.sum(l == 1)), int(np.sum(l < 0))) for l in labs4],
    zref=ZACC, exclus=EXCL, tpr=args.tpr and os.path.relpath(args.tpr, ROOT),
    temps_fin_ps=float(times_ps[-1]), temps_debut_ps=float(times_ps[0]),
)
json.dump(summary, open(os.path.join(OUT, "resume.json"), "w"), indent=1)
