"""
lfpy.py -- interface Python (ctypes) de la bibliotheque C++ libleaflet.so
et fonctions communes aux scripts d'evaluation (lecture des references,
calcul du taux d'accord, construction des fichiers d'index).
"""
import ctypes
import os
import re
import warnings

import numpy as np

warnings.filterwarnings("ignore")
import MDAnalysis as mda  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "Projet3")
if not os.path.isdir(DATA) and os.path.exists(os.path.join(ROOT, "Projet3.zip")):
    # le dossier Projet3/ a ete retire du projet : on lit les donnees de l'archive Projet3.zip,
    # extraite une fois dans un cache hors du dossier du projet
    import zipfile
    DATA = os.path.join(os.path.expanduser("~"), ".cache", "projet_membranes", "Projet3")
    if not os.path.exists(os.path.join(DATA, "large_plasma_membrane", "leaflets.txt")):
        os.makedirs(DATA, exist_ok=True)
        zipfile.ZipFile(os.path.join(ROOT, "Projet3.zip")).extractall(DATA)
RESULTS = os.path.join(ROOT, "resultats")
SYSTEMS = [
    "bilayer_chol",
    "bilayer_peptide",
    "bilayer_prot",
    "model_vesicle",
    "dppc_vesicle",
    "small_plasma_membrane",
    "large_plasma_membrane",
]
VERSIONS = {
    1: "V1 graphe + seuil auto (atome de ref. par nom)",
    2: "V2 graphe + seuil auto (atome de ref. detecte)",
    3: "V3 orientation + normale locale (ACP)",
    4: "V4 finale (orientation lissee)",
    5: "V5 ablation (orientation brute)",
}

_lib = ctypes.CDLL(os.path.join(ROOT, "lib", "libleaflet.so"))
_P = ctypes.c_void_p
_lib.lf_new.restype = _P
_lib.lf_new.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int,
                        ctypes.c_void_p, ctypes.c_char_p]
_lib.lf_free.argtypes = [_P]
_lib.lf_nlip.argtypes = [_P]
_lib.lf_nlip.restype = ctypes.c_int
_lib.lf_molids.argtypes = [_P, ctypes.c_void_p]
_lib.lf_head_atoms.argtypes = [_P, ctypes.c_int, ctypes.c_void_p]
_lib.lf_prepare.argtypes = [_P, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int]
_lib.lf_prepare.restype = ctypes.c_double
_lib.lf_run.argtypes = [_P, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int,
                        ctypes.c_void_p, ctypes.c_void_p]
_lib.lf_run.restype = ctypes.c_int


def _fixed(strings, w=8):
    return b"".join(s.encode()[:w].ljust(w, b" ") for s in strings)


class Analyzer:
    """Topologie (calculee une fois) + analyse d'une trame par le code C++."""

    def __init__(self, universe, nthreads=1, molnums=None, exclude=None):
        """molnums : numero de molecule de chaque atome (p. ex. lu dans un TPR), sinon decoupage par residu.
        exclude : liste de noms de residus a ignorer (p. ex. ["S14"] pour la surfactine)."""
        a = universe.atoms
        self.u = universe
        self.nthreads = nthreads
        resid = np.ascontiguousarray(a.resids, dtype=np.int32)
        self._mol = None if molnums is None else np.ascontiguousarray(molnums, dtype=np.int32)
        ex = ",".join(exclude).encode() if exclude else None
        self._h = _lib.lf_new(len(a), resid.ctypes.data, _fixed(a.resnames), _fixed(a.names), 8,
                              None if self._mol is None else self._mol.ctypes.data, ex)
        self.nlip = _lib.lf_nlip(self._h)
        self.molids = np.zeros(self.nlip, np.int32)
        _lib.lf_molids(self._h, self.molids.ctypes.data)
        self.labels = np.zeros(self.nlip, np.int32)
        self.info = np.zeros(5)
        self.prepared = False

    def __del__(self):
        try:
            _lib.lf_free(self._h)
        except Exception:
            pass

    @staticmethod
    def frame_arrays(ts):
        xyz = np.ascontiguousarray(ts.positions / 10.0, dtype=np.float32)  # A -> nm
        box = np.ascontiguousarray(ts.triclinic_dimensions / 10.0, dtype=np.float32)
        return xyz, box

    def prepare(self, xyz=None, box=None):
        if xyz is None:
            xyz, box = self.frame_arrays(self.u.trajectory.ts)
        t = _lib.lf_prepare(self._h, xyz.ctypes.data, box.ctypes.data, self.nthreads)
        self.prepared = True
        return t

    def head_atoms(self, auto=False):
        out = np.zeros(self.nlip, np.int32)
        _lib.lf_head_atoms(self._h, 1 if auto else 0, out.ctypes.data)
        return out

    def run(self, version=4, cutoff=0.0, xyz=None, box=None):
        if xyz is None:
            xyz, box = self.frame_arrays(self.u.trajectory.ts)
        if version >= 2 and not self.prepared:
            self.prepare(xyz, box)
        _lib.lf_run(self._h, xyz.ctypes.data, box.ctypes.data, version, cutoff, self.nthreads,
                    self.labels.ctypes.data, self.info.ctypes.data)
        return self.labels.copy(), dict(cutoff=self.info[0], time=self.info[1], nleaflets=int(self.info[2]),
                                        unassigned=int(self.info[3]), orphans=int(self.info[4]))


def read_groups(fn):
    """Fichier au format leaflets.txt / .ndx -> liste de (nom, tableau d'entiers)."""
    groups = []
    for line in open(fn):
        m = re.match(r"\s*\[\s*(.+?)\s*\]", line)
        if m:
            groups.append((m.group(1), []))
        elif groups:
            groups[-1][1].extend(int(x) for x in line.split())
    return [(n, np.array(v, int)) for n, v in groups]


def reference_labels(system, molids):
    """Etiquette de reference (0,1,... ou -1 si absent) pour chaque lipide."""
    ref = read_groups(os.path.join(DATA, system, "leaflets.txt"))
    pos = {m: i for i, m in enumerate(molids)}
    lab = np.full(len(molids), -1)
    for k, (_, ids) in enumerate(ref):
        for x in ids:
            if x in pos:
                lab[pos[x]] = k
    return lab


def accuracy(lab, ref):
    """Fraction des lipides de reference places dans le bon feuillet.
    Les feuillets trouves sont associes aux feuillets de reference de facon gloutonne
    (le plus grand d'abord) ; un lipide non assigne (-1) compte comme une erreur."""
    m = ref >= 0
    nref = int(ref.max()) + 1
    found = [L for L in np.unique(lab[lab >= 0])]
    found.sort(key=lambda L: -np.sum(lab == L))
    used, ok = set(), 0
    for L in found:
        c = np.bincount(ref[m & (lab == L)], minlength=nref)
        for k in np.argsort(-c):
            if k not in used:
                used.add(k)
                ok += c[k]
                break
        if len(used) == nref:
            break
    return ok / m.sum()


def write_leaflets(fn, molids, lab):
    with open(fn, "w") as f:
        for k in sorted(set(lab.tolist())):
            ids = molids[lab == k]
            f.write("[ %s ]\n" % ("non_assignes" if k < 0 else "membrane_1_leaflet_%d" % (k + 1)))
            for i in range(0, len(ids), 15):
                f.write("".join("%6d" % x for x in ids[i:i + 15]) + "\n")


def load(system, traj=None):
    gro = os.path.join(DATA, system, system + ".gro")
    return mda.Universe(gro, traj) if traj else mda.Universe(gro)


def tpr_molnums(tpr):
    """Numero de molecule de chaque atome, lu dans un fichier TPR de GROMACS."""
    return mda.Universe(tpr).atoms.molnums


def local_midplane_distance(H, box, cell_nm=3.0, iters=5):
    """Distance signee (nm) de chaque tete au plan median LOCAL d'une bicouche (quasi) plane.
    H : positions des tetes (nm), box : [Lx, Ly, Lz, ...] (nm).
    Le plan median est estime dans chaque case d'une grille de ~cell_nm x cell_nm comme le MILIEU entre la couche
    du haut et celle du bas (separation iterative en deux groupes). Une simple mediane des tetes ne convient pas :
    pour une distribution a deux couches, elle tombe sur l'une des couches des que leurs effectifs different."""
    H = np.asarray(H, float)
    Lx, Ly, Lz = box[0], box[1], box[2]
    z = H[:, 2].copy()
    zc = np.median(z)
    z = zc + (z - zc) - Lz * np.round((z - zc) / Lz)
    # plan median global (deux couches) pour rendre la membrane continue a travers le bord de la boite
    m = z.mean()
    for _ in range(iters):
        up, lo = z[z > m], z[z <= m]
        if len(up) == 0 or len(lo) == 0:
            break
        m = 0.5 * (up.mean() + lo.mean())
    z = m + (z - m) - Lz * np.round((z - m) / Lz)
    nx, ny = max(1, int(Lx / cell_nm)), max(1, int(Ly / cell_nm))
    cell = ((np.floor(H[:, 0] / Lx * nx)).astype(int) % nx) * ny + ((np.floor(H[:, 1] / Ly * ny)).astype(int) % ny)
    zm = np.full(nx * ny, m)
    for c in np.unique(cell):
        zz = z[cell == c]
        mc = m
        for _ in range(iters):
            up, lo = zz[zz > mc], zz[zz <= mc]
            if len(up) == 0 or len(lo) == 0:
                break
            mc = 0.5 * (up.mean() + lo.mean())
        zm[c] = mc
    return z - zm[cell]
