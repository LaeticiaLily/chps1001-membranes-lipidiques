"""
make_tables.py -- tableaux LaTeX du rapport generes a partir des resultats
(resultats/structures/resultats.csv et resultats/trajectoire/resume.json).
Usage : .venv/bin/python scripts/make_tables.py
"""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lfpy import RESULTS, SYSTEMS  # noqa: E402

# le rapport se trouve dans doc_Projet/rapport (ancien emplacement : rapport/)
_ROOT = os.path.dirname(RESULTS)
_REP = os.path.join(_ROOT, "doc_Projet", "rapport")
if not os.path.isdir(_REP):
    _REP = os.path.join(_ROOT, "rapport")
OUT = os.path.join(_REP, "tables")
os.makedirs(OUT, exist_ok=True)
R = list(csv.DictReader(open(os.path.join(RESULTS, "structures", "resultats.csv"))))
NAMES = {"bilayer_chol": "bilayer\\_chol", "bilayer_peptide": "bilayer\\_peptide", "bilayer_prot": "bilayer\\_prot",
         "model_vesicle": "model\\_vesicle", "dppc_vesicle": "dppc\\_vesicle", "small_plasma_membrane": "small\\_PM",
         "large_plasma_membrane": "large\\_PM"}
METHODS = [("V1", "1", "V1"), ("V2", "1", "V2"), ("V3", "1", "V3"), ("V4", "1", "V4"), ("V5", "1", "V5"),
           ("LeafletFinder_15A", "1", "LF 15\\,\\AA"), ("LeafletFinder_optimise", "1", "LF opt."),
           ("FATSLiM", "1", "FATSLiM")]


def get(s, m, t="1"):
    for r in R:
        if r["systeme"] == s and r["methode"] == m and r["threads"] == t:
            return r
    return None


def fr(x, d):
    return ("%.*f" % (d, x)).replace(".", ",")


# --- systemes
with open(os.path.join(OUT, "systemes.tex"), "w") as f:
    for s in SYSTEMS:
        r = get(s, "V4")
        f.write("%s & %s & %s & %s \\\\\n" % (NAMES[s], "{:,}".format(int(r["n_atomes"])).replace(",", "\\,"),
                                               "{:,}".format(int(r["n_lipides"])).replace(",", "\\,"), r["n_ref"]))

# --- accord
with open(os.path.join(OUT, "accord.tex"), "w") as f:
    for s in SYSTEMS:
        cells = []
        vals = []
        for m, t, _ in METHODS:
            r = get(s, m, t)
            vals.append(float(r["accord"]) if r else None)
        best = max(v for v in vals if v is not None)
        for v, (m, t, _) in zip(vals, METHODS):
            r = get(s, m, t)
            if v is None:
                cells.append("échec")
            else:
                c = fr(100 * v, 2)
                if r["methode"].startswith("Leaflet") and int(r["n_feuillets"]) > 2:
                    c += "$^{%s}$" % r["n_feuillets"]
                cells.append("\\textbf{%s}" % c if abs(v - best) < 1e-12 else c)
        f.write("%s & %s \\\\\n" % (NAMES[s], " & ".join(cells)))

# --- erreurs
with open(os.path.join(OUT, "erreurs.tex"), "w") as f:
    for s in SYSTEMS:
        cells = []
        for m, t, _ in METHODS:
            r = get(s, m, t)
            cells.append(r["erreurs"] if r else "--")
        f.write("%s & %s \\\\\n" % (NAMES[s], " & ".join(cells)))

# --- temps
TM = [("V1", "1"), ("V2", "1"), ("V3", "1"), ("V4", "1"), ("V4", "12"), ("LeafletFinder_15A", "1"),
      ("LeafletFinder_optimise", "1"), ("FATSLiM", "1"), ("FATSLiM", "12")]
with open(os.path.join(OUT, "temps.tex"), "w") as f:
    for s in SYSTEMS:
        cells = []
        for m, t in TM:
            r = get(s, m, t)
            if not r:
                cells.append("échec")
                continue
            v = float(r["t_min_ms"])
            d = 3 if v < 1 else (2 if v < 10 else (1 if v < 100 else 0))
            cells.append(fr(v, d))
        f.write("%s & %s \\\\\n" % (NAMES[s], " & ".join(cells)))

# --- seuils automatiques
with open(os.path.join(OUT, "seuils.tex"), "w") as f:
    for s in SYSTEMS:
        r1, r2, ro = get(s, "V1"), get(s, "V2"), get(s, "LeafletFinder_optimise")
        f.write("%s & %s & %s & %s & %s \\\\\n" % (
            NAMES[s], fr(float(r1["seuil_nm"]), 2), fr(float(r2["seuil_nm"]), 2),
            fr(float(ro["seuil_nm"]), 2) if ro else "échec", fr(float(get(s, "V2")["t_prep_ms"]), 1)))

# --- trajectoires
TRAJS = [("trajectoire_surfactine_20us", "Surf. prod\\_dt1000ns.xtc"),
         ("trajectoire_surfactine_prod1", "Surf. prod1.xtc"),
         ("trajectoire_surfactine_prodtrr", "Surf. prod.trr (tronqué)"),
         ("trajectoire", "Vésicule DPPC (test)")]
Js = []
for d, name in TRAJS:
    fn = os.path.join(RESULTS, d, "resume.json")
    if os.path.exists(fn):
        Js.append((d, name, json.load(open(fn))))


def dur(J):
    t = (J.get("temps_fin_ps", J["duree_ps"]) - J.get("temps_debut_ps", 0.0))
    return ("%s\\,µs" % fr(t / 1e6, 1)) if t >= 1e5 else ("%s\\,ns" % fr(t / 1e3, 1))


with open(os.path.join(OUT, "traj_donnees.tex"), "w") as f:
    for d, name, J in Js:
        f.write("%s & %d & %s & %s & %s \\\\\n" % (name, J["n_trames"], dur(J), J["n_lipides"],
                                                   "{:,}".format(J["n_atomes"]).replace(",", "\\,")))
with open(os.path.join(OUT, "traj_temps.tex"), "w") as f:
    for d, name, J in Js:
        n = J["n_trames"]
        f.write("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n" % (
            name, fr(J["V4_total_ms"], 0 if J["V4_total_ms"] > 100 else 1), fr(J["V4_ms_trame"], 2),
            fr(J["V1_ms_trame"], 2), fr(J["LF_ms_trame"], 1), fr(J["FATSLiM1_total_ms"] / n, 1),
            fr(J["debit_par_ms_trame"], 3), fr(J["extrapol_100k_par_s"], 0)))
with open(os.path.join(OUT, "traj_accord.tex"), "w") as f:
    for d, name, J in Js:
        z = J.get("zref") or {}
        zz = lambda k: fr(100 * z[k][0], 2) if k in z else "--"
        f.write("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n" % (
            name, fr(100 * J["accord_V4_FATSLiM_moy"], 2), fr(100 * J["accord_V4_V1_moy"], 2),
            fr(100 * J["accord_V4_LF_moy"], 2), zz("V4"), zz("V1"), zz("FATSLiM"), zz("LeafletFinder")))
with open(os.path.join(OUT, "trajectoire_macros.tex"), "w") as f:
    def m(k, v):
        f.write("\\newcommand{\\%s}{%s}\n" % (k, v))
    for tag, d in (("P", "trajectoire_surfactine_prod1"), ("L", "trajectoire_surfactine_20us"), ("D", "trajectoire")):
        fn = os.path.join(RESULTS, d, "resume.json")
        if not os.path.exists(fn):
            continue
        J = json.load(open(fn))
        m("T%sTrames" % tag, J["n_trames"])
        m("T%sLipides" % tag, J["n_lipides"])
        m("T%sAtomes" % tag, "{:,}".format(J["n_atomes"]).replace(",", "\\,"))
        m("T%sLecture" % tag, fr(J["lecture_ms"] / J["n_trames"], 2))
        m("T%sDetect" % tag, fr(J["detection_tetes_ms"], 1))
        m("T%sVq" % tag, fr(J["V4_ms_trame"], 2))
        m("T%sVtot" % tag, fr(J["V4_total_ms"] / 1000.0, 2))
        m("T%sSeq" % tag, fr(J["debit_seq_ms_trame"], 3))
        m("T%sPar" % tag, fr(J["debit_par_ms_trame"], 3))
        m("T%sSeqHz" % tag, "%d" % round(1000 / J["debit_seq_ms_trame"]))
        m("T%sParHz" % tag, "%d" % round(1000 / J["debit_par_ms_trame"]))
        m("T%sExtSeq" % tag, fr(J["extrapol_100k_seq_s"], 0))
        m("T%sExtPar" % tag, fr(J["extrapol_100k_par_s"], 0))
        m("T%sFat" % tag, fr(J["FATSLiM1_total_ms"] / 1000.0, 1))
        m("T%sLF" % tag, fr(J["LF_total_ms"] / 1000.0, 1))
        m("T%sAccFat" % tag, fr(100 * J["accord_V4_FATSLiM_moy"], 2))
        m("T%sAccFatMin" % tag, fr(100 * J["accord_V4_FATSLiM_min"], 2))
        m("T%sChgMoy" % tag, "%d" % round(np.mean(J["changements"])) if J["changements"] else "0")
        m("T%sAccRef" % tag, fr(100 * J["accord_ref_trame0"], 2) if J.get("accord_ref_trame0") is not None else "n/a")
        nf2 = sum(1 for t in J["tailles"] if min(t[0], t[1]) < 0.25 * (t[0] + t[1]))
        m("T%sFusions" % tag, nf2)
    m("TNcpu", os.cpu_count())
# style du rapport (comme le compte-rendu du TP7) : un filet horizontal sous chaque ligne
for fn in os.listdir(OUT):
    if fn.endswith(".tex") and "macros" not in fn:
        t = open(os.path.join(OUT, fn)).read()
        t = "\n".join((l + " \\hline") if l.rstrip().endswith("\\\\") else l for l in t.split("\n"))
        open(os.path.join(OUT, fn), "w").write(t)
print("tableaux ecrits dans", OUT)
