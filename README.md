# Projet CHPS1001 — Analyse de la structure des membranes lipidiques

Groupe : BERNARD Valentin, BRAHAM Amyne, FOURNAISE Antoine, ZAMOUM Tassadit — Master 2 CHPS, URCA, 2026-2027.

Programme qui range chaque lipide d'une membrane simulée dans son feuillet (couche du haut ou du bas,
intérieur ou extérieur d'une vésicule), écrit en C++/OpenMP, comparé aux outils existants LeafletFinder
(MDAnalysis) et FATSLiM, et appliqué à une trajectoire.

## Par où commencer (rien à installer pour lire)

| Ordre | Fichier | Ce que c'est |
|---|---|---|
| 1 | [`rapport/explications.pdf`](rapport/explications.pdf) | **Explications simples** de tout le projet, point par point, avec les fichiers correspondants |
| 2 | [`doc_Projet/rapport/rapport.pdf`](doc_Projet/rapport/rapport.pdf) | **Le rapport rendu** (21 pages) |
| 3 | [`resultats/verification/resume_verification.txt`](resultats/verification/resume_verification.txt) | Vérification lipide par lipide : pour chaque membrane et chaque méthode, lipides bien placés / mal placés |
| 4 | [`resultats/figures/`](resultats/figures) et [`fig_verification.png`](resultats/verification/fig_verification.png) | Les figures |

**Résultats en bref** : la version finale (V4) place correctement 99,77 à 100 % des lipides des 7 membranes
(les rares erreurs sont des cholestérols au milieu de la membrane), aussi bien que FATSLiM mais 28 à 250 fois
plus vite. Sur la trajectoire du cours, elle analyse une trame de 3040 lipides en 1,3 ms.

## Qui est quoi dans le dépôt

| Dossier / fichier | Contenu |
|---|---|
| `src/leaflet.cpp` | **Tout le calcul** : les versions V1 à V5 (option `-m 1` à `-m 5`) |
| `Makefile` | Compilation (`make`) → `bin/leaflet` (programme) et `lib/libleaflet.so` (utilisé par Python) |
| `scripts/benchmark_structures.py` | Points 1 à 4 et 6 : nos versions, LeafletFinder et FATSLiM sur les 7 membranes |
| `scripts/run_trajectory.py` | Points 5 et 6 : analyse d'une trajectoire, image par image, et temps de calcul |
| `scripts/verifier_resultats.py` | Comparaison lipide par lipide avec les bonnes réponses (`leaflets.txt`) |
| `scripts/make_figures.py`, `make_tables.py` | Figures et tableaux du rapport |
| `scripts/lfpy.py` | Lien entre Python et le programme C++ ; calcul du taux de réussite |
| `Projet3/` | Les 7 membranes fournies (`.gro`) et leurs bonnes réponses (`leaflets.txt`) |
| `trajectoire/` | Petite trajectoire de test (vésicule DPPC, 11 images) |
| `resultats/structures/` | Résultats sur les 7 membranes (`resultats.csv`, feuillets trouvés par chaque méthode) |
| `resultats/trajectoire_surfactine_*/` | Résultats sur la trajectoire du cours (`resume.txt` = bilan lisible) |
| `resultats/verification/` | Vérification lipide par lipide (CSV, figure, scripts VMD) |
| `doc_Projet/rapport/` | Le rapport (PDF + sources LaTeX) ; `rapport/` : les explications |
| `installer.sh` | Installe tout ce qu'il faut (voir ci-dessous) |
| `resultats/relance_tout.sh` | Relance tous les calculs, la vérification, les figures et les PDF |

## Installer et relancer (Linux, sans droits administrateur)

```bash
git clone https://github.com/LaeticiaLily/chps1001-membranes-lipidiques.git
cd chps1001-membranes-lipidiques
bash installer.sh                 # ~5 min : Python, MDAnalysis, FATSLiM, compilation du C++
bash resultats/relance_tout.sh    # ~15 min : tous les calculs, la vérification, les figures et les PDF
```

La trajectoire du cours (`Surfactine-R1/`, fournie par l'enseignant) n'est pas publiée dans ce dépôt :
pour refaire son analyse, copier le dossier `Surfactine-R1` du cours à la racine du projet. Sans lui, le script
de relance saute simplement ces étapes ; les résultats déjà calculés sont dans `resultats/`.

## Vérifier soi-même

```bash
./bin/leaflet -f Projet3/bilayer_chol/bilayer_chol.gro -m 4 -r Projet3/bilayer_chol/leaflets.txt
.venv/bin/python scripts/verifier_resultats.py
vmd -e resultats/verification/vmd_large_plasma_membrane.tcl     # lipides mal placés en rouge
```

## Remarques

- Les temps ont été mesurés sur un portable sur batterie : `bash resultats/relance_tout.sh` sur secteur,
  applications fermées, donne des temps plus nets (le rapport se met à jour tout seul).
- `Surfactine-R1/prod.trr` (fourni) est tronqué : 61 images lisibles. `prod1.xtc` (1001 images) sert de mesure principale.
