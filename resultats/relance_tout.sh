#!/bin/bash
# =============================================================================
#  Relance COMPLETE du projet : compilation, calculs, verification, figures, PDF
#  Usage (depuis n'importe ou) :  bash "resultats/relance_tout.sh"
#  Conseil : brancher le portable sur secteur et fermer les autres applications
#  (les temps mesures en dependent). Duree : environ 15 minutes.
# =============================================================================
set -e
cd "$(dirname "$0")/.."
PY=.venv/bin/python
S="--gro Surfactine-R1/prod.gro --tpr Surfactine-R1/prod.tpr --exclude S14 --zref"
etape() { echo; echo "=================== $1 ==================="; }

etape "1. Compilation du programme C++ (bin/leaflet, lib/libleaflet.so)"
make

etape "2. Les 7 structures : V1..V5, LeafletFinder, FATSLiM (points 1 a 4 et 6)"
$PY scripts/benchmark_structures.py

etape "3. Trajectoire de test (vesicule DPPC, 11 trames)"
$PY scripts/run_trajectory.py --out trajectoire

# une trajectoire absente (fichiers volumineux non copies) est sautee, sans arreter le reste
si() { if [ -f "$1" ]; then return 0; else echo "  -> $1 absent : etape sautee"; return 1; fi; }
etape "4. Trajectoire du cours : production 20 us (21 trames)"
si Surfactine-R1/prod_dt1000ns.xtc && $PY scripts/run_trajectory.py $S --traj Surfactine-R1/prod_dt1000ns.xtc --out trajectoire_surfactine_20us
etape "5. Trajectoire du cours : prod1.xtc (1001 trames)"
si Surfactine-R1/prod1.xtc && $PY scripts/run_trajectory.py $S --traj Surfactine-R1/prod1.xtc --out trajectoire_surfactine_prod1
etape "6. Trajectoire du cours : prod.trr (61 trames lisibles, fichier tronque)"
si Surfactine-R1/prod.trr && $PY scripts/run_trajectory.py $S --traj Surfactine-R1/prod.trr --frames 0:61 --out trajectoire_surfactine_prodtrr

etape "7. Changements de feuillet par type de lipide (trajectoire du cours)"
si Surfactine-R1/prod_dt1000ns.xtc && $PY scripts/analyse_types_trajectoire.py Surfactine-R1/prod.gro \
    Surfactine-R1/prod_dt1000ns.xtc Surfactine-R1/prod.tpr resultats/trajectoire_surfactine_20us/analyse_types.txt S14
si Surfactine-R1/prod1.xtc && $PY scripts/analyse_types_trajectoire.py Surfactine-R1/prod.gro \
    Surfactine-R1/prod1.xtc Surfactine-R1/prod.tpr resultats/trajectoire_surfactine_prod1/analyse_types.txt S14

etape "8. VERIFICATION lipide par lipide (resultats/verification/)"
$PY scripts/verifier_resultats.py

etape "9. Figures et tableaux du rapport"
$PY scripts/make_figures.py
$PY scripts/make_tables.py

etape "10. Rapport et explications en PDF"
(cd doc_Projet/rapport && for i in 1 2; do pdflatex -interaction=nonstopmode rapport.tex >/dev/null; done; pdflatex -interaction=nonstopmode rapport.tex | grep "Output written")
(cd rapport && for i in 1 2; do pdflatex -interaction=nonstopmode explications.tex >/dev/null; done; pdflatex -interaction=nonstopmode explications.tex | grep "Output written")
rm -f doc_Projet/rapport/rapport.{aux,log,out,toc} rapport/explications.{aux,log,out,toc}

etape "TERMINE"
echo "Resume de la verification : resultats/verification/resume_verification.txt"
