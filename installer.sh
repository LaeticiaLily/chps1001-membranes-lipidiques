#!/bin/bash
# =============================================================================
#  Installation de l'environnement du projet (sans droits administrateur)
#  Usage, depuis le dossier du projet :  bash installer.sh
#
#  Prerequis deja presents sur la plupart des Linux : g++ (avec OpenMP), make,
#  python3 (>= 3.10), pdflatex (pour recompiler les PDF, facultatif).
#  Cree : .venv/ (MDAnalysis, SciPy, Matplotlib...) et .venv-fatslim/ (FATSLiM,
#  qui ne fonctionne qu'avec Python 3.9 et NumPy < 1.24), puis compile le C++.
# =============================================================================
set -e
cd "$(dirname "$0")"

echo "=== 1. Environnement Python principal (.venv)"
python3 -m venv --system-site-packages .venv
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q numpy scipy matplotlib networkx MDAnalysis

echo "=== 2. FATSLiM dans un environnement Python 3.9 (.venv-fatslim)"
if ! command -v uv >/dev/null 2>&1; then
    echo "  installation de l'outil uv (gestionnaire de versions de Python)"
    python3 -m pip install --user -q uv || .venv/bin/pip install -q uv
    export PATH="$HOME/.local/bin:$PWD/.venv/bin:$PATH"
fi
uv python install 3.9
uv venv -q --python 3.9 .venv-fatslim
uv pip install -q --python .venv-fatslim/bin/python "numpy<1.24" "cython<3" "setuptools<60" wheel pip
.venv-fatslim/bin/python -m pip install -q --no-build-isolation fatslim

echo "=== 3. Compilation du programme C++"
make

echo
echo "Installation terminee. Tests rapides :"
./bin/leaflet -f Projet3/bilayer_chol/bilayer_chol.gro -m 4 -r Projet3/bilayer_chol/leaflets.txt | grep "Accord"
.venv-fatslim/bin/python .venv-fatslim/bin/fatslim --version | head -3
echo "Pour tout relancer : bash resultats/relance_tout.sh"
