"""
make_figures.py -- figures du rapport a partir des CSV de resultats/ (style : scripts/figstyle.py).
Usage : .venv/bin/python scripts/make_figures.py
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import C, INK, INK2, RED, WIDTH, origin, plt, save  # noqa: E402
from lfpy import RESULTS, SYSTEMS, Analyzer, load, local_midplane_distance  # noqa: E402

FIG = os.path.join(RESULTS, "figures")
os.makedirs(FIG, exist_ok=True)
SHORT = {"bilayer_chol": "bilayer_chol", "bilayer_peptide": "bilayer_peptide", "bilayer_prot": "bilayer_prot",
         "model_vesicle": "model_vesicle", "dppc_vesicle": "dppc_vesicle", "small_plasma_membrane": "small_PM",
         "large_plasma_membrane": "large_PM"}


def rows(fn):
    p = os.path.join(RESULTS, fn)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else None


def fl(r, k):
    try:
        return float(r[k])
    except (ValueError, KeyError, TypeError):
        return np.nan


R = rows("structures/resultats.csv")


def get(system, method, threads="1", key="t_min_ms"):
    for r in R:
        if r["systeme"] == system and r["methode"] == method and r["threads"] == threads:
            return fl(r, key)
    return np.nan


def centered(P, box, center_on=None):
    """Positions centrees sur la vesicule (image minimale), pour que l'origine soit son centre."""
    L = box[:3]
    sel = slice(None) if center_on is None else center_on
    d0 = P - P[0]
    d0 -= L * np.round(d0 / L)
    c = P[0] + d0[sel].mean(0)
    v = P - c
    return v - L * np.round(v / L)


# ================================================================ 1. balayage du seuil (V1)
scan = rows("structures/balayage_seuil.csv")
fig, axes = plt.subplots(4, 2, figsize=(WIDTH, 11.0))
for ax, s in zip(axes.flat, SYSTEMS):
    d = [r for r in scan if r["systeme"] == s]
    x = np.array([fl(r, "seuil_nm") for r in d])
    ax.plot(x, [fl(r, "taille_1") for r in d], color=C[0], label="plus grand groupe")
    ax.plot(x, [fl(r, "taille_2") for r in d], color=C[1], label="deuxième groupe")
    ax.plot(x, [fl(r, "accord") for r in d], color=INK2, lw=1.4, ls="--", label="accord avec la référence")
    cut = get(s, "V1", key="seuil_nm")
    ax.axvline(cut, color=C[2], lw=2)
    ax.text(cut + 0.06, 0.16, ("%.2f nm" % cut).replace(".", ","), color=INK, fontsize=10, ha="left", va="bottom",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.5, alpha=0.9))
    ax.set_title(SHORT[s])
    ax.set_xlim(0, 2.6)
    ax.set_ylim(0, 1.05)
    origin(ax)
for ax in list(axes.flat)[:7]:
    ax.set_xlabel("seuil (nm)")
    ax.set_ylabel("fraction des lipides")
axes.flat[-1].axis("off")
h, l = axes.flat[0].get_legend_handles_labels()
axes.flat[-1].legend(h + [plt.Line2D([], [], color=C[2], lw=2)], l + ["seuil choisi par V1"], loc="center")
fig.tight_layout()
save(fig, os.path.join(FIG, "fig_balayage_seuil.png"))

# ================================================================ 2. temps d'analyse des structures
methods = [("V1", "V1"), ("V3", "V3"), ("V4", "V4 (finale)"), ("FATSLiM", "FATSLiM"),
           ("LeafletFinder_15A", "LeafletFinder 15 Å")]
fig, ax = plt.subplots(figsize=(WIDTH, 4.4))
w = 0.16
x = np.arange(len(SYSTEMS))
for k, (m, lab) in enumerate(methods):
    ax.bar(x + (k - 2) * w, [get(s, m) for s in SYSTEMS], w * 0.92, color=C[k], label=lab)
ax.set_yscale("log")
ax.set_ylim(0.01, 3000)
ax.set_xlim(-0.5, len(SYSTEMS) - 0.5)
ax.set_xticks(x, [SHORT[s] for s in SYSTEMS], rotation=25, ha="right")
ax.set_ylabel("temps par structure (ms, log)")
ax.grid(axis="x", visible=False)
ax.legend(ncol=3, loc="upper left", bbox_to_anchor=(0, 1.22))
fig.tight_layout()
save(fig, os.path.join(FIG, "fig_temps_structures.png"))

# ================================================================ 3. scalabilite OpenMP
sc = rows("structures/scalabilite_threads.csv")
fig, ax = plt.subplots(figsize=(6.0, 3.8))
for v, col in (("1", C[0]), ("4", C[2])):
    d = [r for r in sc if r["version"] == v]
    ax.plot([int(r["threads"]) for r in d], [fl(r, "t_min_ms") for r in d], color=col, marker="o", ms=6,
            label="V%s" % v)
ax.set_xlabel("nombre de threads OpenMP")
ax.set_ylabel("temps (ms)")
ax.set_title("large_plasma_membrane (19 239 lipides)")
ax.set_xlim(0, 12.5)
ax.set_xticks(range(0, 13, 2))
origin(ax)
ax.legend(loc="lower left")
fig.tight_layout()
save(fig, os.path.join(FIG, "fig_threads.png"))

# ================================================================ 4. trajectoires : temps par trame
for d, out in (("trajectoire_surfactine_prod1", "fig_trajectoire_temps.png"),
               ("trajectoire", "fig_trajectoire_dppc_temps.png")):
    T = rows(os.path.join(d, "par_trame.csv"))
    if not T:
        continue
    fig, ax = plt.subplots(figsize=(WIDTH, 4.0))
    tr = [int(r["trame"]) for r in T]
    many = len(T) > 60
    for col, lab, c in (("LeafletFinder_ms", "LeafletFinder 15 Å", C[4]), ("FATSLiM_ms_1thr", "FATSLiM", C[3]),
                        ("V1_ms", "V1", C[0]), ("V4_ms", "V4 (finale)", C[2])):
        ax.plot(tr, [fl(r, col) for r in T], color=c, lw=1.3 if many else 2.2, marker=None if many else "o",
                ms=5, label=lab)
    ax.set_yscale("log")
    ax.set_ylim(0.5, 300)
    ax.set_xlim(0, tr[-1])
    origin(ax, ylog=True)
    ax.set_xlabel("trame")
    ax.set_ylabel("temps par trame (ms, log)")
    ax.legend(ncol=4, loc="upper left", bbox_to_anchor=(0, 1.16))
    fig.tight_layout()
    save(fig, os.path.join(FIG, out))

# ================================================================ 5. production de 20 us
T = rows("trajectoire_surfactine_20us/par_trame.csv")
if T:
    t = np.array([fl(r, "temps_ps") for r in T]) / 1e6
    fig, axes = plt.subplots(3, 1, figsize=(WIDTH, 9.6))
    ax = axes[0]
    diff = np.array([fl(r, "V4_feuillet1") - fl(r, "V4_feuillet2") for r in T])
    ax.bar(t, diff, width=0.6, color=C[0])
    ax.set_ylabel("feuillet 1 − feuillet 2\n(lipides)")
    ax.set_title("(a) Écart de taille entre les deux feuillets (≈ 1520 lipides chacun)", pad=12)
    lim = max(5, np.abs(diff).max() * 1.2)
    ax.set_ylim(-lim, lim)
    ax.set_xlim(0, 20.6)
    ax.set_xlabel("temps (µs)")
    origin(ax, ysigned=True)
    ax = axes[1]
    ax.bar(t[1:], [fl(r, "changements_V4_vs_trame_prec") for r in T[1:]], width=0.6, color=C[2])
    ax.set_ylabel("lipides ayant\nchangé de feuillet")
    ax.set_title("(b) Changements de feuillet depuis la trame précédente :\nuniquement du cholestérol")
    ax.set_xlabel("temps (µs)")
    ax.set_xlim(0, 20.6)
    origin(ax)
    ax = axes[2]
    for col, lab, c, ls in (("zref_V4", "V4", C[2], "-"), ("zref_V1", "V1", C[0], "--"),
                            ("zref_FATSLiM", "FATSLiM", C[3], ":"), ("zref_LeafletFinder", "LeafletFinder", C[4], "-")):
        ax.plot(t, [100 * fl(r, col) for r in T], color=c, ls=ls, marker="o", ms=4, label=lab)
    ax.set_ylabel("accord (%)")
    ax.set_xlabel("temps (µs)")
    ax.set_title("(c) Accord avec la référence géométrique\n(tête au-dessus ou au-dessous du plan médian)")
    ax.set_xlim(0, 20.6)
    ax.set_ylim(0, 105)
    origin(ax)
    ax.legend(ncol=4, loc="lower right")
    fig.tight_layout()
    save(fig, os.path.join(FIG, "fig_surfactine_20us.png"))

# ================================================================ 6. vues des feuillets trouves par V4
fig = plt.figure(figsize=(WIDTH, 8.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1])
for k, s in enumerate(["dppc_vesicle", "model_vesicle"]):
    ax = fig.add_subplot(gs[0, k])
    u = load(s)
    A = Analyzer(u)
    lab, _ = A.run(4)
    hp = u.atoms.positions[A.head_atoms()] / 10.0
    box = u.dimensions / 10.0
    P = centered(hp, box, lab >= 0) if s == "dppc_vesicle" else hp - hp.mean(0)
    m = np.abs(P[:, 2]) < 1.5
    for L, col, name in ((0, C[0], "feuillet 1 (extérieur)"), (1, C[1], "feuillet 2 (intérieur)"),
                         (-1, INK2, "non assigné")):
        sel = m & (lab == L)
        ax.scatter(P[sel, 0], P[sel, 1], s=7, color=col, lw=0, label=name)
    lim = np.abs(P[m & (lab >= 0), :2]).max() * 1.12  # cadre limite a la membrane (pas aux lipides libres)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_title("%s : coupe |z| < 1,5 nm" % s)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.grid(False)
    origin(ax, xsigned=True, ysigned=True)
ax = fig.add_subplot(gs[1, :])
u = load("large_plasma_membrane")
A = Analyzer(u)
lab, _ = A.run(4)
hp = u.atoms.positions[A.head_atoms()] / 10.0
d = local_midplane_distance(hp, u.dimensions[:3] / 10.0)
m = hp[:, 1] < 8.0
rn = u.atoms.resnames[A.head_atoms()]
for L, col in ((0, C[0]), (1, C[1])):
    sel = m & (lab == L) & (rn != "CHOL")
    ax.scatter(hp[sel, 0], d[sel], s=5, color=col, lw=0)
    sel = m & (lab == L) & (rn == "CHOL")
    ax.scatter(hp[sel, 0], d[sel], s=16, color=col, marker="^", lw=0)
ax.set_xlim(0, u.dimensions[0] / 10.0)
ax.set_ylim(-3.2, 3.2)
ax.set_xlabel("x (nm)")
ax.set_ylabel("distance au plan\nmédian (nm)")
ax.set_title("large_plasma_membrane : vue de côté (tranche y < 8 nm) ; triangles = cholestérol")
ax.grid(False)
origin(ax, ysigned=True)
h, l = fig.axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, markerscale=3, bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.03, 1, 1))
save(fig, os.path.join(FIG, "fig_feuillets.png"))

# ================================================================ 7. membrane de la trajectoire surfactine
SURF = os.path.join(os.path.dirname(RESULTS), "Surfactine-R1")
if os.path.exists(os.path.join(SURF, "prod.tpr")):
    import MDAnalysis as mda
    from lfpy import tpr_molnums
    u = mda.Universe(os.path.join(SURF, "prod.gro"), os.path.join(SURF, "prod_dt1000ns.xtc"))
    u.trajectory[-1]
    A = Analyzer(u, molnums=tpr_molnums(os.path.join(SURF, "prod.tpr")), exclude=["S14"])
    lab, _ = A.run(4)
    h = A.head_atoms(auto=True)
    P = u.atoms.positions[h] / 10.0
    box = u.dimensions[:3] / 10.0
    d = local_midplane_distance(P, box)
    zmid = np.median(P[:, 2] - d)  # plan median moyen, pour placer les surfactines
    rn = u.atoms.resnames[h]
    surf = u.select_atoms("resname S14 and name SC3").positions / 10.0
    fig, ax = plt.subplots(figsize=(WIDTH, 4.2))
    m = P[:, 1] < 6.0
    for k, (col, name) in enumerate([(C[0], "feuillet 1"), (C[1], "feuillet 2")]):
        sel = m & (lab == k) & (rn != "CHOL")
        ax.scatter(P[sel, 0], d[sel], s=12, color=col, lw=0, label=name + " : PIPC, GCER")
        sel = m & (lab == k) & (rn == "CHOL")
        ax.scatter(P[sel, 0], d[sel], s=26, color=col, marker="^", lw=0, label=name + " : cholestérol")
    ms = surf[:, 1] < 6.0
    ax.scatter(surf[ms, 0], surf[ms, 2] - zmid, s=40, color=INK, marker="x", lw=1.4, label="surfactine (queue)",
               zorder=6)
    ax.set_xlim(0, box[0])
    ax.set_ylim(-4.5, 4.5)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("distance au plan\nmédian (nm)")
    ax.set_title("Trajectoire surfactine, t = 20 µs : vue de côté (tranche y < 6 nm)")
    ax.grid(False)
    origin(ax, ysigned=True)
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2), markerscale=1.3)
    fig.tight_layout()
    save(fig, os.path.join(FIG, "fig_surfactine_vue.png"))
print("figures ecrites dans", FIG)
