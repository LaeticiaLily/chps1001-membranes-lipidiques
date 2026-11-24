# Projet CHPS1001 — Analyse de la structure des membranes lipidiques

Groupe : BERNARD Valentin, BRAHAM Amyne, FOURNAISE Antoine, ZAMOUM Tassadit — Master 2 CHPS, URCA, 2026-2027.

Assignation des lipides d'une membrane à leur feuillet (C++/OpenMP), comparaison avec
LeafletFinder (MDAnalysis) et FATSLiM.

**Rapport : `doc_Projet/rapport/rapport.pdf`** (source `rapport.tex` au même endroit).
Explications détaillées pour la réunion : `rapport/explications.pdf`.

| Élément | Emplacement |
|---|---|
| Énoncé, données | `2026_CHPS1001_...pdf`, `Projet3/` |
| Code C++ (versions 1 à 5) | `src/leaflet.cpp`, `Makefile` → `bin/leaflet`, `lib/libleaflet.so` |
| Scripts Python | `scripts/` (évaluation, trajectoire, figures, tableaux) |
| Résultats | `resultats/structures/`, `resultats/trajectoire/`, `resultats/figures/` |
| Trajectoire du cours | `Surfactine-R1/` (GROMACS, Martini 3 : 3040 lipides + 96 surfactines) |
| Trajectoire de test | `trajectoire/` (vésicule DPPC, 11 trames, issue des données de test de FATSLiM) |
| Environnements Python | `.venv/` (MDAnalysis…), `.venv-fatslim/` (FATSLiM 0.2.2, Python 3.9) |

## Installation (Linux, sans droits administrateur)
```bash
git clone https://github.com/LaeticiaLily/chps1001-membranes-lipidiques.git
cd chps1001-membranes-lipidiques
bash installer.sh            # environnements Python (MDAnalysis, FATSLiM) + compilation du C++
```
La trajectoire du cours (`Surfactine-R1/`, fournie par l'enseignant) n'est pas publiée dans ce dépôt :
copier le dossier du cours à la racine du projet pour refaire l'analyse de la trajectoire.
Les 7 structures (`Projet3/`) et la trajectoire de test (`trajectoire/`) sont incluses.

## Tout relancer
```bash
cd ~/Documents/"projet membranne lipidique"
bash resultats/relance_tout.sh          # compilation, calculs, vérification, figures, PDF (~15 min, sur secteur)
```

## Vérifier les résultats lipide par lipide
```bash
.venv/bin/python scripts/verifier_resultats.py             # ~1 min, une fois les calculs faits
less resultats/verification/resume_verification.txt        # bons / mauvais / non assignés, par méthode
vmd -e resultats/verification/vmd_large_plasma_membrane.tcl  # erreurs de V4 en rouge dans VMD
./bin/leaflet -f Projet3/bilayer_chol/bilayer_chol.gro -m 4 -r Projet3/bilayer_chol/leaflets.txt -o mes_feuillets.txt
```
Fichiers : `resultats/verification/<membrane>_V4_lipides.csv` (une ligne par lipide : attendu / trouvé),
`<membrane>_<méthode>_erreurs.csv` (lipides mal placés), `fig_verification.png`.

## À faire
- Facultatif : relancer `resultats/relance_tout.sh` sur secteur, applications fermées, pour des temps plus nets
  (les mesures actuelles ont été faites sur batterie), puis `scripts/make_figures.py`, `scripts/make_tables.py`
  et recompiler le rapport.
- `Surfactine-R1/prod.trr` est tronqué (61 trames lisibles) : s'il existe une copie complète, la remettre à la place.
