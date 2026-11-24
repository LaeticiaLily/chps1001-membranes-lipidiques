"""Convertit une trajectoire (xtc/trr) en fichier GRO multi-trames lisible par bin/leaflet.
(le GROWriter de MDAnalysis n'ecrit qu'une trame par fichier : les trames sont concatenees)
Usage : .venv/bin/python scripts/xtc_to_gro.py structure.gro trajectoire.xtc sortie.gro"""
import os
import sys
import tempfile
import warnings
warnings.filterwarnings("ignore")
import MDAnalysis as mda
u = mda.Universe(sys.argv[1], sys.argv[2])
with open(sys.argv[3], "w") as out, tempfile.TemporaryDirectory() as tmp:
    fn = os.path.join(tmp, "f.gro")
    for ts in u.trajectory:
        u.atoms.write(fn, reindex=False)
        out.write(open(fn).read())
print(len(u.trajectory), "trames ->", sys.argv[3])
