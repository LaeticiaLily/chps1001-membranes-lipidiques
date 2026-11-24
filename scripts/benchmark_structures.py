"""
benchmark_structures.py -- evaluation des versions V1..V5 et comparaison avec
LeafletFinder (MDAnalysis) et FATSLiM sur les 7 structures fournies.

Sorties (dossier resultats/structures) :
  resultats.csv                         une ligne par (systeme, methode)
  <systeme>/<methode>_leaflets.txt      feuillets trouves (format leaflets.txt)
  balayage_seuil.csv                    taille des 2 plus grandes composantes vs seuil
  scalabilite_threads.csv               temps de V4 selon le nombre de threads

Usage : .venv/bin/python scripts/benchmark_structures.py [--quick]
"""
import csv
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lfpy import (DATA, RESULTS, ROOT, SYSTEMS, VERSIONS, Analyzer, accuracy, load,  # noqa: E402
                  reference_labels, write_leaflets)
from MDAnalysis.analysis.leaflet import LeafletFinder, optimize_cutoff  # noqa: E402
from MDAnalysis.lib.distances import self_capped_distance  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

QUICK = "--quick" in sys.argv
OUT = os.path.join(RESULTS, "structures")
FATSLIM = [os.path.join(ROOT, ".venv-fatslim", "bin", "python"), os.path.join(ROOT, ".venv-fatslim", "bin", "fatslim")]
NCPU = os.cpu_count()


def timed_runs(fn, min_reps=5, min_time=0.5, max_reps=200):
    """Execute fn plusieurs fois ; renvoie (dernier resultat, temps min, temps median) en secondes."""
    times, res = [], None
    t_start = time.perf_counter()
    while len(times) < min_reps or (time.perf_counter() - t_start < min_time and len(times) < max_reps):
        t0 = time.perf_counter()
        res = fn()
        times.append(time.perf_counter() - t0)
    return res, min(times), float(np.median(times))


def labels_from_groups(groups, nlip):
    """Liste de tableaux d'indices de lipides -> etiquettes (groupes tries par taille)."""
    lab = np.full(nlip, -1)
    for k, g in enumerate(sorted(groups, key=len, reverse=True)):
        lab[g] = k
    return lab


def run_leafletfinder(u, heads, head2lip, cutoff_A, optimize=False):
    ag = u.atoms[heads]

    def go():
        c = cutoff_A
        if optimize:
            c, _ = optimize_cutoff(u, ag, dmin=8.0, dmax=22.0, step=0.5, max_imbalance=0.2, pbc=True)
        lf = LeafletFinder(u, ag, cutoff=c, pbc=True)
        return lf, c

    reps = 1 if optimize else 3
    (lf, c), tmin, tmed = timed_runs(go, min_reps=reps, min_time=0.0, max_reps=reps)
    groups = [head2lip[g.ix] for g in lf.groups()]
    return labels_from_groups(groups, int(head2lip.max()) + 1), tmin, tmed, c / 10.0, len(groups)


def run_fatslim(gro, heads, head2lip_dict, nthreads, nlip, traj=None):
    gro = os.path.abspath(gro)  # FATSLiM est lance dans un dossier temporaire
    traj = os.path.abspath(traj) if traj else None
    with tempfile.TemporaryDirectory() as tmp:
        ndx = os.path.join(tmp, "hg.ndx")
        h = heads + 1
        with open(ndx, "w") as f:
            f.write("[ headgroups ]\n")
            for i in range(0, len(h), 15):
                f.write(" ".join(map(str, h[i:i + 15])) + "\n")
        cmd = FATSLIM + ["membranes", "-c", gro, "-n", ndx, "--output-index-hg", os.path.join(tmp, "hg_out.ndx"),
                         "--nthreads", str(nthreads)]
        if traj:
            # FATSLiM nomme ses sorties d'apres le temps de la trame (arrondi) : des trames rapprochees
            # s'ecraseraient. On lui donne une copie de la trajectoire ou la trame k est au temps k ns.
            import MDAnalysis as mda
            u_ = mda.Universe(gro, traj)
            tx = os.path.join(tmp, "traj_renum.xtc")
            with mda.Writer(tx, n_atoms=len(u_.atoms)) as w:
                for k, ts in enumerate(u_.trajectory):
                    ts.time = 1000.0 * k
                    w.write(u_.atoms)
            cmd += ["-t", tx]
        t0 = time.perf_counter()
        p = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
        wall = time.perf_counter() - t0
        frame_ms = [float(x) for x in re.findall(r"done in (\d+) ms", p.stdout)]
        # un fichier par trame, nomme d'apres le temps (hg_out_0000.ndx, hg_out_1000.ndx...) : tri numerique
        outs = sorted((f for f in os.listdir(tmp) if f.startswith("hg_out")),
                      key=lambda f: float(re.findall(r"_(\d+)\.ndx$", f)[0]) if re.findall(r"_(\d+)\.ndx$", f) else 0)
        labs = []
        for fn in outs:
            lab = np.full(nlip, -1)
            k = -1
            for line in open(os.path.join(tmp, fn)):
                if "[" in line:
                    k += 1
                    continue
                for x in line.split():
                    lab[head2lip_dict[int(x) - 1]] = k
            labs.append(lab)
        return labs, frame_ms, wall, p.stdout + p.stderr


def threshold_scan(P, box, ref, name, rows):
    """Taille des deux plus grandes composantes connexes selon le seuil (graphe des atomes de reference)."""
    n = len(P)
    rmax = min(2.6, 0.49 * min(box[:3]) / 10.0 * 10 / 10)  # nm
    pairs, d = self_capped_distance(P * 10.0, rmax * 10.0, box=box)
    d = d / 10.0
    for c in np.arange(0.6, rmax + 1e-9, 0.02):
        m = d < c
        A = coo_matrix((np.ones(m.sum()), (pairs[m, 0], pairs[m, 1])), shape=(n, n))
        _, lab = connected_components(A, directed=False)
        sizes = np.sort(np.bincount(lab))[::-1]
        s2 = sizes[1] if len(sizes) > 1 else 0
        rows.append(dict(systeme=name, seuil_nm=round(c, 3), n_composantes=len(sizes), taille_1=sizes[0] / n,
                         taille_2=s2 / n, accord=accuracy(lab, ref)))


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, scan_rows = [], []
    systems = SYSTEMS[:3] if QUICK else SYSTEMS
    for s in systems:
        print("=" * 70, "\n", s, flush=True)
        os.makedirs(os.path.join(OUT, s), exist_ok=True)
        gro = os.path.join(DATA, s, s + ".gro")
        u = load(s)
        A = Analyzer(u, nthreads=1)
        ref = reference_labels(s, A.molids)
        nref = int((ref >= 0).sum())
        xyz, box = A.frame_arrays(u.trajectory.ts)
        t_prep = A.prepare(xyz, box)
        base = dict(systeme=s, n_lipides=A.nlip, n_ref=nref, n_atomes=len(u.atoms))

        # --- nos versions
        for v in VERSIONS:
            (lab, info), tmin, tmed = timed_runs(lambda: A.run(v, xyz=xyz, box=box))
            acc = accuracy(lab, ref)
            write_leaflets(os.path.join(OUT, s, "V%d_leaflets.txt" % v), A.molids, lab)
            sizes = sorted(np.bincount(lab[lab >= 0]).tolist(), reverse=True)
            rows.append(dict(base, methode="V%d" % v, threads=1, accord=acc, erreurs=int(round((1 - acc) * nref)),
                             n_feuillets=info["nleaflets"], tailles=" ".join(map(str, sizes[:4])),
                             non_assignes=int((lab < 0).sum()), seuil_nm=round(info["cutoff"], 3),
                             t_min_ms=tmin * 1e3, t_med_ms=tmed * 1e3,
                             t_prep_ms=t_prep * 1e3 if v >= 2 else 0.0))
            print("  V%d  acc=%.4f  t=%.3f ms  seuil=%.2f  feuillets=%s" % (v, acc, tmin * 1e3, info["cutoff"], sizes[:3]),
                  flush=True)
        # V4 multi-threads
        A12 = Analyzer(u, nthreads=NCPU)
        A12.prepare(xyz, box)
        (lab, info), tmin, tmed = timed_runs(lambda: A12.run(4, xyz=xyz, box=box))
        rows.append(dict(base, methode="V4", threads=NCPU, accord=accuracy(lab, ref), n_feuillets=info["nleaflets"],
                         non_assignes=int((lab < 0).sum()), seuil_nm=info["cutoff"], t_min_ms=tmin * 1e3,
                         t_med_ms=tmed * 1e3, erreurs=int(round((1 - accuracy(lab, ref)) * nref))))

        heads = A.head_atoms(auto=False)
        head2lip = np.full(len(u.atoms), -1)
        head2lip[heads] = np.arange(A.nlip)
        hd = {int(h): i for i, h in enumerate(heads)}

        # --- balayage du seuil (figure)
        threshold_scan(xyz[heads], u.dimensions, ref, s, scan_rows)

        # --- LeafletFinder (MDAnalysis)
        for label, cut, opt in [("LeafletFinder_15A", 15.0, False), ("LeafletFinder_optimise", None, True)]:
            try:
                lab, tmin, tmed, c, ng = run_leafletfinder(u, heads, head2lip, cut, optimize=opt)
                acc = accuracy(lab, ref)
                sizes = sorted(np.bincount(lab[lab >= 0]).tolist(), reverse=True)
                write_leaflets(os.path.join(OUT, s, "%s_leaflets.txt" % label), A.molids, lab)
                rows.append(dict(base, methode=label, threads=1, accord=acc, erreurs=int(round((1 - acc) * nref)),
                                 n_feuillets=ng, tailles=" ".join(map(str, sizes[:4])), non_assignes=int((lab < 0).sum()),
                                 seuil_nm=c, t_min_ms=tmin * 1e3, t_med_ms=tmed * 1e3))
                print("  %-24s acc=%.4f  t=%.1f ms  seuil=%.2f  groupes=%d %s" % (label, acc, tmin * 1e3, c, ng, sizes[:3]),
                      flush=True)
            except Exception as e:  # pragma: no cover
                print("  %s : ECHEC %s" % (label, e))

        # --- FATSLiM
        for nt in (1, NCPU):
            labs, fms, wall, log = run_fatslim(gro, heads, hd, nt, A.nlip)
            if not labs:
                print("  FATSLiM ECHEC\n", log[-500:])
                continue
            lab = labs[0]
            acc = accuracy(lab, ref)
            sizes = sorted(np.bincount(lab[lab >= 0]).tolist(), reverse=True)
            if nt == 1:
                write_leaflets(os.path.join(OUT, s, "FATSLiM_leaflets.txt"), A.molids, lab)
            rows.append(dict(base, methode="FATSLiM", threads=nt, accord=acc, erreurs=int(round((1 - acc) * nref)),
                             n_feuillets=int(lab.max() + 1), tailles=" ".join(map(str, sizes[:4])),
                             non_assignes=int((lab < 0).sum()), seuil_nm=2.0, t_min_ms=fms[0] if fms else np.nan,
                             t_med_ms=fms[0] if fms else np.nan, t_total_processus_ms=wall * 1e3))
            print("  FATSLiM (%2d thr)          acc=%.4f  t=%s ms (processus %.0f ms) feuillets=%s" %
                  (nt, acc, fms, wall * 1e3, sizes[:3]), flush=True)

    keys = ["systeme", "n_atomes", "n_lipides", "n_ref", "methode", "threads", "accord", "erreurs", "n_feuillets",
            "tailles", "non_assignes", "seuil_nm", "t_min_ms", "t_med_ms", "t_prep_ms", "t_total_processus_ms"]
    with open(os.path.join(OUT, "resultats.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (("%.6g" % r[k]) if isinstance(r.get(k), float) else r.get(k, "")) for k in keys})
    with open(os.path.join(OUT, "balayage_seuil.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(scan_rows[0].keys()))
        w.writeheader()
        w.writerows(scan_rows)

    # --- scalabilite OpenMP de V4 sur la plus grande membrane
    if not QUICK:
        s = "large_plasma_membrane"
        u = load(s)
        with open(os.path.join(OUT, "scalabilite_threads.csv"), "w") as f:
            f.write("threads,version,t_min_ms,t_med_ms\n")
            for nt in [1, 2, 4, 6, 8, 10, 12]:
                A = Analyzer(u, nthreads=nt)
                xyz, box = A.frame_arrays(u.trajectory.ts)
                A.prepare(xyz, box)
                for v in (1, 4):
                    _, tmin, tmed = timed_runs(lambda: A.run(v, xyz=xyz, box=box), min_reps=20, min_time=1.0)
                    f.write("%d,%d,%.4f,%.4f\n" % (nt, v, tmin * 1e3, tmed * 1e3))
                    print("threads=%2d V%d t=%.3f ms" % (nt, v, tmin * 1e3), flush=True)
    print("Resultats ecrits dans", OUT)


if __name__ == "__main__":
    main()
