"""
verifier_resultats.py -- verification LIPIDE PAR LIPIDE des feuillets trouves, par rapport aux
fichiers de reference fournis (Projet3/<systeme>/leaflets.txt).

Pour chaque systeme et chaque methode (V1..V5, LeafletFinder, FATSLiM), on relit le fichier de feuillets
produit par scripts/benchmark_structures.py (resultats/structures/<systeme>/<methode>_leaflets.txt)
et on le compare a la reference. On recalcule aussi V4 a l'instant, pour verifier que le fichier enregistre
correspond bien au programme actuel.

Sorties (resultats/verification/) :
  resume_verification.txt        tableau : bons / mauvais feuillet / non assignes / groupes parasites, par methode
  <systeme>_V4_lipides.csv       une ligne PAR LIPIDE : reference, trouve, statut, position, orientation
  <systeme>_<methode>_erreurs.csv  les seuls lipides mal places, pour chaque methode
  fig_verification.png           ou sont les lipides mal places par V4 (croix rouges)
  vmd_<systeme>.tcl              visualisation dans VMD : vmd -e resultats/verification/vmd_<systeme>.tcl

Usage : .venv/bin/python scripts/verifier_resultats.py
"""
import collections
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lfpy import DATA, RESULTS, ROOT, SYSTEMS, Analyzer, load, read_groups  # noqa: E402
from MDAnalysis.lib.distances import minimize_vectors  # noqa: E402

from figstyle import C as PAL, RED, WIDTH, origin, plt  # noqa: E402

OUT = os.path.join(RESULTS, "verification")
os.makedirs(OUT, exist_ok=True)
METHODS = ["V1", "V2", "V3", "V4", "V5", "LeafletFinder_15A", "LeafletFinder_optimise", "FATSLiM"]
VESICLES = ("model_vesicle", "dppc_vesicle")


def groups_to_labels(groups, molids):
    """Fichier de feuillets -> etiquette par lipide (-1 = non assigne). Les groupes sont tries par taille."""
    pos = {m: i for i, m in enumerate(molids)}
    lab = np.full(len(molids), -1)
    gs = [(n, ids) for n, ids in groups if not n.startswith("non_assign")]
    gs.sort(key=lambda g: -len(g[1]))
    for k, (_, ids) in enumerate(gs):
        for x in ids:
            if x in pos:
                lab[pos[x]] = k
    return lab, [n for n, _ in gs]


def match(lab, ref, nref):
    """Associe chaque groupe trouve a un feuillet de reference (le plus grand groupe d'abord).
    Renvoie un tableau 'feuillet de reference attribue' par lipide (-1 : non assigne, -2 : groupe parasite)."""
    out = np.full(len(lab), -1)
    used = set()
    groups = sorted(set(lab[lab >= 0].tolist()), key=lambda g: -np.sum(lab == g))
    for g in groups:
        c = np.bincount(ref[(lab == g) & (ref >= 0)], minlength=nref)
        k = next((k for k in np.argsort(-c) if k not in used), None)
        if k is None:
            out[lab == g] = -2
        else:
            used.add(k)
            out[lab == g] = k
    return out


def molecule_atoms(u):
    """Atomes de chaque molecule (suite d'atomes de meme residu), comme dans le programme C++."""
    a = u.atoms
    ch = np.r_[True, (a.resids[1:] != a.resids[:-1]) | (a.resnames[1:] != a.resnames[:-1])]
    st = np.flatnonzero(ch)
    en = np.r_[st[1:], len(a)]
    return st, en


def geometry(u, A, system):
    """Position de la tete, distance a la surface mediane et orientation, pour chaque lipide."""
    heads = A.head_atoms(auto=False)
    X = u.atoms.positions / 10.0
    box = u.dimensions.copy()
    box[:3] /= 10.0
    st, en = molecule_atoms(u)
    start_of = {int(s): (int(s), int(e)) for s, e in zip(st, en)}
    # le debut de la molecule de chaque lipide = plus grand debut <= atome de tete
    idx = np.searchsorted(st, heads, side="right") - 1
    H = X[heads]
    o = np.zeros_like(H)
    for i, k in enumerate(idx):
        s, e = start_of[int(st[k])]
        d = minimize_vectors(X[s:e] - H[i], box)
        c = d.mean(0)
        n = np.linalg.norm(c)
        o[i] = -c / n if n > 0 else 0
    return heads, H, o, box


def surface_distance(system, H, ref, box):
    """Distance signee de la tete a la surface mediane de la membrane, et direction 'vers l'exterieur'."""
    if system in VESICLES:
        L = box[:3]
        d0 = H - H[0]
        d0 -= L * np.round(d0 / L)
        c = H[0] + d0[ref >= 0].mean(0)
        v = H - c
        v -= L * np.round(v / L)
        r = np.linalg.norm(v, axis=1)
        rmid = 0.5 * (np.median(r[ref == 0]) + np.median(r[ref == 1]))
        return r - rmid, v / r[:, None], c
    # bicouche plane : plan median local (milieu entre les deux couches de tetes, grille d'environ 3 nm)
    from lfpy import local_midplane_distance
    nrm = np.zeros_like(H)
    nrm[:, 2] = 1.0
    return local_midplane_distance(H, box[:3]), nrm, None


STATUT = {0: "correct", 1: "mauvais_feuillet", 2: "non_assigne", 3: "groupe_parasite", 4: "absent_de_la_reference"}


def statuses(found_ref, ref):
    s = np.zeros(len(ref), int)
    s[(ref >= 0) & (found_ref >= 0) & (found_ref != ref)] = 1
    s[(ref >= 0) & (found_ref == -1)] = 2
    s[(ref >= 0) & (found_ref == -2)] = 3
    s[ref < 0] = 4
    return s


def write_vmd(system, gro, molids, ref, st4):
    """Script VMD : feuillet 1 de reference en bleu, feuillet 2 en orange, erreurs de V4 en rouge (grosses spheres)."""
    def sel(ids):
        ids = [int(x) - 1 for x in ids]  # numero de molecule (1, 2...) -> 'residue' VMD (0, 1...)
        return "residue " + " ".join(map(str, ids)) if ids else "none"
    fn = os.path.join(OUT, "vmd_%s.tcl" % system)
    with open(fn, "w") as f:
        f.write("# Verification des feuillets de %s dans VMD\n" % system)
        f.write("# Bleu / orange : feuillets de REFERENCE (leaflets.txt) ; rouge : lipides que V4 place AUTREMENT.\n")
        f.write("# Lancer depuis le dossier du projet : vmd -e %s\n" % os.path.relpath(fn, ROOT))
        f.write("mol new {%s} waitfor all\n" % os.path.relpath(gro, ROOT))
        f.write("mol delrep 0 top\n")
        for k, col in ((0, "blue"), (1, "orange")):
            f.write("mol representation Points 3\nmol color ColorID %d\nmol selection {%s}\nmol addrep top\n"
                    % (0 if col == "blue" else 3, sel(molids[ref == k])))
        bad = molids[(st4 == 1) | (st4 == 2) | (st4 == 3)]
        f.write("mol representation VDW 1.4 12\nmol color ColorID 1\nmol selection {%s}\nmol addrep top\n" % sel(bad))
        f.write("display projection Orthographic\naxes location off\ncatch {package require pbctools; pbc box}\n")
        f.write("puts \"%s : %d lipides mal places par V4 (en rouge)\"\n" % (system, len(bad)))
    return fn


def main():
    summary = []
    fig, axes = plt.subplots(4, 2, figsize=(WIDTH, 12.5))
    for ax, system in zip(axes.flat, SYSTEMS):
        gro = os.path.join(DATA, system, system + ".gro")
        u = load(system)
        A = Analyzer(u)
        A.prepare()
        molids = A.molids
        refg = read_groups(os.path.join(DATA, system, "leaflets.txt"))
        ref, refnames = groups_to_labels(refg, molids)
        # l'ordre de reference doit rester celui du fichier, pas celui de la taille
        pos = {m: i for i, m in enumerate(molids)}
        ref = np.full(len(molids), -1)
        for k, (_, ids) in enumerate(refg):
            for x in ids:
                ref[pos[x]] = k
        refnames = [n for n, _ in refg]
        nref = len(refg)
        heads, H, o, box = geometry(u, A, system)
        dist, nrm, _ = surface_distance(system, H, ref, box)
        cosn = (o * nrm).sum(1)  # orientation du lipide par rapport a la normale (+1 : tete vers l'exterieur/haut)
        rn = u.atoms.resnames[heads]

        # V4 recalcule maintenant
        lab_now, _ = A.run(4)
        res_now = match(lab_now, ref, nref)
        for meth in METHODS:
            fn = os.path.join(RESULTS, "structures", system, "%s_leaflets.txt" % meth)
            if not os.path.exists(fn):
                continue
            lab, _ = groups_to_labels(read_groups(fn), molids)
            fref = match(lab, ref, nref)
            st = statuses(fref, ref)
            m = ref >= 0
            cnt = collections.Counter(st[m])
            bytype = collections.Counter(rn[m & (st != 0)])
            summary.append(dict(systeme=system, methode=meth, n_ref=int(m.sum()), corrects=cnt[0],
                                mauvais_feuillet=cnt[1], non_assignes=cnt[2], groupes_parasites=cnt[3],
                                hors_reference=int((ref < 0).sum()), accord=cnt[0] / m.sum(),
                                erreurs_par_type=" ".join("%s:%d" % kv for kv in bytype.most_common())))
            with open(os.path.join(OUT, "%s_%s_erreurs.csv" % (system, meth)), "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["numero_molecule", "type", "feuillet_reference", "feuillet_trouve", "statut",
                            "x_nm", "y_nm", "z_nm", "distance_surface_mediane_nm", "orientation_cos"])
                for i in np.flatnonzero(m & (st != 0)):
                    w.writerow([molids[i], rn[i], refnames[ref[i]],
                                refnames[fref[i]] if fref[i] >= 0 else STATUT[st[i]], STATUT[st[i]],
                                "%.2f" % H[i, 0], "%.2f" % H[i, 1], "%.2f" % H[i, 2], "%+.2f" % dist[i], "%+.2f" % cosn[i]])
            if meth == "V4":
                st4 = st
                same = np.array_equal(fref, res_now)
                summary[-1]["V4_fichier_identique_au_calcul_actuel"] = "oui" if same else "NON"
                with open(os.path.join(OUT, "%s_V4_lipides.csv" % system), "w", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(["numero_molecule", "type", "feuillet_reference", "feuillet_trouve_V4", "statut",
                                "x_nm", "y_nm", "z_nm", "distance_surface_mediane_nm", "orientation_cos"])
                    for i in range(len(molids)):
                        w.writerow([molids[i], rn[i], refnames[ref[i]] if ref[i] >= 0 else "absent",
                                    refnames[fref[i]] if fref[i] >= 0 else STATUT[min(st[i], 3) if st[i] else 2],
                                    STATUT[st[i]], "%.2f" % H[i, 0], "%.2f" % H[i, 1], "%.2f" % H[i, 2],
                                    "%+.2f" % dist[i], "%+.2f" % cosn[i]])
        write_vmd(system, gro, molids, ref, st4)

        # figure : ou sont les erreurs de V4 ? (origine a l'intersection des axes)
        bad = (ref >= 0) & (st4 != 0)
        if system in VESICLES:
            L = box[:3]
            d0 = H - H[0]
            d0 -= L * np.round(d0 / L)
            c = H[0] + d0[ref >= 0].mean(0)
            P = (H - c) - L * np.round((H - c) / L) if system != "model_vesicle" else H - H.mean(0)
            sl = np.abs(P[:, 2]) < 1.5
            for k, col in ((0, PAL[0]), (1, PAL[1])):
                s_ = sl & (ref == k)
                ax.scatter(P[s_, 0], P[s_, 1], s=4, color=col, lw=0)
            ax.scatter(P[bad, 0], P[bad, 1], s=70, marker="x", color=RED, lw=2, zorder=6)
            lim = np.abs(P[sl & (ref >= 0), :2]).max() * 1.15
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_aspect("equal")
            ax.set_xlabel("x (nm)")
            ax.set_ylabel("y (nm)")
            ax.set_title("%s : %d erreur(s) / %d" % (system, bad.sum(), (ref >= 0).sum()))
            origin(ax, xsigned=True, ysigned=True)
        else:
            for k, col in ((0, PAL[0]), (1, PAL[1])):
                s_ = ref == k
                ax.scatter(H[s_, 0], dist[s_], s=3 if len(H) > 5000 else 6, color=col, lw=0, alpha=0.7)
            ax.scatter(H[bad, 0], dist[bad], s=70, marker="x", color=RED, lw=2, zorder=6)
            ax.set_xlim(0, box[0])
            ax.set_ylim(-3.5, 3.5)
            ax.set_xlabel("x (nm)")
            ax.set_ylabel("distance au plan\nmédian (nm)")
            ax.set_title("%s : %d erreur(s) / %d" % (system.replace("_plasma_membrane", "_PM"), bad.sum(),
                                                      (ref >= 0).sum()))
            origin(ax, ysigned=True)
        print("%-22s V4 : %d erreurs sur %d (fichier enregistre identique au calcul actuel : %s)" %
              (system, bad.sum(), (ref >= 0).sum(), summary[[s["methode"] for s in summary].index("V4", len(summary) - 8)]
               .get("V4_fichier_identique_au_calcul_actuel")), flush=True)
    ax = axes.flat[-1]
    ax.axis("off")
    ax.scatter([], [], color=PAL[0], s=40, label="feuillet 1 de référence")
    ax.scatter([], [], color=PAL[1], s=40, label="feuillet 2 de référence")
    ax.scatter([], [], color=RED, marker="x", s=70, lw=2, label="lipide placé autrement\npar V4")
    ax.legend(loc="center", frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_verification.png"))

    with open(os.path.join(OUT, "resume_verification.txt"), "w") as f:
        f.write("VERIFICATION LIPIDE PAR LIPIDE PAR RAPPORT AUX FICHIERS leaflets.txt FOURNIS\n")
        f.write("corrects = dans le bon feuillet ; mauvais = dans l'autre feuillet ; non_ass. = non assigne ;\n")
        f.write("parasites = dans un petit groupe en plus des 2 feuillets ; hors_ref = lipides absents de la reference\n\n")
        cur = None
        for s in summary:
            if s["systeme"] != cur:
                cur = s["systeme"]
                f.write("=== %s (%d lipides de reference)\n" % (cur, s["n_ref"]))
                f.write("  %-24s %9s %8s %9s %10s %8s %9s   %s\n" % ("methode", "corrects", "mauvais", "non_ass.",
                                                                  "parasites", "hors_ref", "accord", "erreurs par type"))
            f.write("  %-24s %9d %8d %9d %10d %8d %8.3f%%   %s\n" % (
                s["methode"], s["corrects"], s["mauvais_feuillet"], s["non_assignes"], s["groupes_parasites"],
                s["hors_reference"], 100 * s["accord"], s["erreurs_par_type"]))
        f.write("\nV4 : le fichier enregistre est identique au calcul refait maintenant pour : %s\n" % ", ".join(
            s["systeme"] for s in summary if s.get("V4_fichier_identique_au_calcul_actuel") == "oui"))
    print(open(os.path.join(OUT, "resume_verification.txt")).read())


def verify_trajectory():
    """Trajectoire du cours (pas de leaflets.txt) : comparaison lipide par lipide, sur la premiere et la derniere trame
    de la production, de V4 avec FATSLiM, V1 et la reference geometrique (plan median local)."""
    import MDAnalysis as mda
    from lfpy import tpr_molnums, local_midplane_distance
    from benchmark_structures import run_fatslim
    surf = os.path.join(ROOT, "Surfactine-R1")
    gro, xtc, tpr = (os.path.join(surf, f) for f in ("prod.gro", "prod_dt1000ns.xtc", "prod.tpr"))
    if not all(os.path.exists(f) for f in (gro, xtc, tpr)):
        return
    u = mda.Universe(gro, xtc)
    A = Analyzer(u, molnums=tpr_molnums(tpr), exclude=["S14"])
    A.prepare()
    hn = A.head_atoms(auto=False)
    ha = A.head_atoms(auto=True)
    rn = u.atoms.resnames[ha]
    labsF, _, _, _ = run_fatslim(gro, hn, {int(x): i for i, x in enumerate(hn)}, os.cpu_count(), A.nlip, traj=xtc)
    lines = ["", "VERIFICATION SUR LA TRAJECTOIRE DU COURS (Surfactine-R1/prod_dt1000ns.xtc, %d lipides)" % A.nlip,
             "Pas de fichier de reference : on compare V4, lipide par lipide, a FATSLiM, a V1 et a la position de la",
             "tete par rapport au plan median local (au-dessus = feuillet du haut)."]
    for k in (0, len(u.trajectory) - 1):
        ts = u.trajectory[k]
        l4, _ = A.run(4)
        l1, _ = A.run(1)
        H = u.atoms.positions[hn] / 10
        d = local_midplane_distance(H, ts.dimensions[:3] / 10)
        zr = (d > 0).astype(int)
        al = lambda x: x if np.mean(x == l4) >= 0.5 else np.where(x >= 0, 1 - x, x)
        lf, l1, zr = al(labsF[k]), al(l1), al(zr)
        lines.append("  trame %d (t = %.0f us) :" % (k, ts.time / 1e6))
        for name, other in (("FATSLiM", lf), ("V1", l1), ("reference geometrique", zr)):
            dif = l4 != other
            lines.append("    V4 different de %-22s : %4d lipides (%.2f %%)  types : %s" % (
                name, dif.sum(), 100 * dif.mean(), dict(collections.Counter(rn[dif]).most_common())))
        fn = os.path.join(OUT, "trajectoire_surfactine_trame%02d_lipides.csv" % k)
        with open(fn, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["numero_molecule", "type", "V4", "FATSLiM", "V1", "reference_geometrique",
                        "distance_plan_median_nm", "V4_identique_aux_3_autres"])
            for i in range(A.nlip):
                w.writerow([A.molids[i], rn[i], l4[i] + 1, lf[i] + 1, l1[i] + 1, zr[i] + 1, "%+.2f" % d[i],
                            "oui" if (l4[i] == lf[i] == l1[i] == zr[i]) else "NON"])
        dd = d[l4 != zr]
        if len(dd):
            lines.append("    lipides ou V4 et la reference geometrique different : |distance au plan median| "
                         "mediane %.2f nm (les feuillets sont a ~2 nm)" % np.median(np.abs(dd)))
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "resume_verification.txt"), "a") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
    verify_trajectory()
