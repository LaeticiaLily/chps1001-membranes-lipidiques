"""Ecrit le decoupage en molecules d'un fichier TPR au format lu par bin/leaflet (--molsizes) :
une ligne "nombre_de_molecules atomes_par_molecule" par bloc de molecules identiques.
Usage : .venv/bin/python scripts/tpr_molsizes.py topol.tpr molsizes.txt"""
import sys
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import MDAnalysis as mda
mn = mda.Universe(sys.argv[1]).atoms.molnums
sizes = np.bincount(mn)
with open(sys.argv[2], "w") as f:
    i = 0
    while i < len(sizes):
        j = i
        while j < len(sizes) and sizes[j] == sizes[i]:
            j += 1
        f.write("%d %d\n" % (j - i, sizes[i]))
        i = j
print(len(sizes), "molecules ->", sys.argv[2])
