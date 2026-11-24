"""
figstyle.py -- style commun des figures (lisibles une fois inserees en pleine largeur dans le rapport A4).

Regle des axes : l'origine (0, 0) est a l'intersection des axes.
  - grandeur positive : l'axe commence exactement a 0 (pas de marge avant 0) ;
  - grandeur signee (distance au plan median, coordonnee centree sur une vesicule) : l'axe des abscisses
    (ou des ordonnees) passe par 0, au milieu du graphe.
Les echelles logarithmiques (temps de calcul) ne contiennent pas 0 : seule l'autre coordonnee part de 0.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# palette categorielle (ordre fixe), encre et grille discretes
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
RED = "#e34948"
WIDTH = 7.5  # pouces : une figure de 7,5" inseree sur 16 cm garde des textes d'environ 9 pt

plt.rcParams.update({
    "font.size": 11, "axes.titlesize": 11.5, "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.7, "axes.axisbelow": True, "legend.frameon": False,
    "savefig.dpi": 220, "savefig.bbox": "tight", "lines.linewidth": 2.2, "axes.linewidth": 1.0,
})


def _comma(x, pos=None):
    """Graduations a la francaise : virgule decimale, pas de zeros inutiles."""
    s = ("%.6g" % x).replace(".", ",").replace("-", "\u2212")
    return s


def origin(ax, xsigned=False, ysigned=False, xlog=False, ylog=False):
    """Place l'origine a l'intersection des axes (voir la regle en tete du fichier).
    A appeler APRES avoir fixe les limites. Pour un axe signe, l'axe passe par 0 au milieu du graphe,
    son nom est ecrit au bout de l'axe (comme un repere mathematique) et les graduations ont un fond blanc
    pour rester lisibles au-dessus des points."""
    from matplotlib.ticker import FuncFormatter
    if not xlog:
        ax.xaxis.set_major_formatter(FuncFormatter(_comma))
        if not xsigned:
            ax.set_xlim(left=0)
    if not ylog:
        ax.yaxis.set_major_formatter(FuncFormatter(_comma))
        if not ysigned:
            ax.set_ylim(bottom=0)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    if xsigned and not xlog:
        ax.spines["left"].set_position(("data", 0))
        fx = (0 - x0) / (x1 - x0)
        lab = ax.get_ylabel()
        ax.set_ylabel(lab, rotation=0, ha="left", va="top")
        ax.yaxis.set_label_coords(fx + 0.03, 1.0)  # au bout de l'axe, juste a droite
    if ysigned and not ylog:
        ax.spines["bottom"].set_position(("data", 0))
        fy = (0 - y0) / (y1 - y0)
        lab = ax.get_xlabel()
        ax.set_xlabel(lab, ha="right", va="top")
        ax.xaxis.set_label_coords(1.0, fy - 0.09)  # au bout de l'axe, juste en dessous
    if xsigned and ysigned:
        # un seul « 0 » a l'origine
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, p: "" if abs(x) < 1e-12 else _comma(x)))
    if xsigned or ysigned:
        ax.set_axisbelow(False)  # axes et graduations au-dessus des points et des barres
        ax.grid(False)
        ax.figure.canvas.draw()
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            t.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.9))


def save(fig, path):
    fig.savefig(path)
    plt.close(fig)
