// =============================================================================
//  leaflet.cpp  --  Assignation des lipides d'une membrane a leur feuillet
//  Projet CHPS1001 (URCA, M2 CHPS) : Analyse de la structure des membranes
//
//  Versions de l'algorithme (option -m) :
//    1 : graphe sur atomes de reference (noms connus) + seuil automatique
//    2 : idem, mais l'atome de reference est detecte automatiquement
//        (atome le plus "expose" de chaque type de lipide)
//    3 : orientation des lipides + normale locale par ACP (type FATSLiM)
//    4 : version finale : orientation lissee sur les voisins (sans ACP),
//        seuil fixe physique, rattachement des orphelins
//    5 : ablation : comme 4 mais orientation brute (sans lissage)
//
//  Compilation : voir Makefile (g++ -O3 -march=native -fopenmp)
//  Le meme fichier fournit l'executable (main) et la bibliotheque C
//  (libleaflet.so, API lf_*) utilisee depuis Python via ctypes.
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <malloc.h>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using std::string;
using std::vector;

// Les tableaux de travail (quelques Mo) sont realloues a chaque trame : on demande a la glibc
// de les garder dans le tas au lieu de les rendre au systeme (mmap/munmap + remise a zero
// des pages par le noyau), ce qui coutait ~30 % du temps d'analyse.
static void tune_malloc() {
    static bool done = false;
    if (done) return;
    done = true;
    mallopt(M_MMAP_THRESHOLD, 1 << 30);
    mallopt(M_TRIM_THRESHOLD, 1 << 30);
}

static inline double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ----------------------------------------------------------------------------
// Boite de simulation (triclinique generale, vecteurs en lignes : r = f . H)
// ----------------------------------------------------------------------------
struct Box {
    double H[3][3], Hi[3][3], w[3];
    bool ortho;
    void set(const double h[3][3]) {
        memcpy(H, h, sizeof(H));
        ortho = (H[0][1] == 0 && H[0][2] == 0 && H[1][0] == 0 && H[1][2] == 0 && H[2][0] == 0 && H[2][1] == 0);
        double det = H[0][0] * (H[1][1] * H[2][2] - H[1][2] * H[2][1]) - H[0][1] * (H[1][0] * H[2][2] - H[1][2] * H[2][0]) +
                     H[0][2] * (H[1][0] * H[2][1] - H[1][1] * H[2][0]);
        Hi[0][0] = (H[1][1] * H[2][2] - H[1][2] * H[2][1]) / det;
        Hi[0][1] = (H[0][2] * H[2][1] - H[0][1] * H[2][2]) / det;
        Hi[0][2] = (H[0][1] * H[1][2] - H[0][2] * H[1][1]) / det;
        Hi[1][0] = (H[1][2] * H[2][0] - H[1][0] * H[2][2]) / det;
        Hi[1][1] = (H[0][0] * H[2][2] - H[0][2] * H[2][0]) / det;
        Hi[1][2] = (H[0][2] * H[1][0] - H[0][0] * H[1][2]) / det;
        Hi[2][0] = (H[1][0] * H[2][1] - H[1][1] * H[2][0]) / det;
        Hi[2][1] = (H[0][1] * H[2][0] - H[0][0] * H[2][1]) / det;
        Hi[2][2] = (H[0][0] * H[1][1] - H[0][1] * H[1][0]) / det;
        for (int k = 0; k < 3; k++)  // largeur perpendiculaire selon l'axe fractionnaire k
            w[k] = 1.0 / std::sqrt(Hi[0][k] * Hi[0][k] + Hi[1][k] * Hi[1][k] + Hi[2][k] * Hi[2][k]);
    }
    inline void frac(const float* r, double* f) const {
        for (int k = 0; k < 3; k++) f[k] = r[0] * Hi[0][k] + r[1] * Hi[1][k] + r[2] * Hi[2][k];
    }
    // image minimale d'un petit vecteur (intra-moleculaire)
    inline void minimg(float* d) const {
        if (ortho) {
            for (int k = 0; k < 3; k++) d[k] -= (float)(H[k][k] * std::nearbyint(d[k] / H[k][k]));
        } else {
            double f[3];
            frac(d, f);
            for (int k = 0; k < 3; k++) f[k] -= std::nearbyint(f[k]);
            for (int k = 0; k < 3; k++) d[k] = (float)(f[0] * H[0][k] + f[1] * H[1][k] + f[2] * H[2][k]);
        }
    }
};

// ----------------------------------------------------------------------------
// Topologie : quelles molecules sont des lipides, de quel type, atomes de tete
// ----------------------------------------------------------------------------
static const char* NONLIPID[] = {"W",   "SOL", "HOH", "WAT", "TIP3", "PW",  "NA+", "CL-", "NA",  "CL",  "K",   "K+",  "CA",
                                 "CA+", "ION", "MG",  "ZN",  "ACE",  "NH2", "NME", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN",
                                 "GLU", "GLY", "HIS", "HSD", "HSE",  "HIE", "HID", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
                                 "SER", "THR", "TRP", "TYR", "VAL",  nullptr};
// noms usuels de l'atome de reference (phosphate ou equivalent), par priorite
static const char* HEADNAMES[] = {"PO4", "P8", "P", "CP", "ROH", "O3", "AM1", "GL1", nullptr};

struct Topo {
    int natoms = 0, nlip = 0;
    vector<int> start, len, type, molid;  // par lipide (molid : rang 1..N de la molecule dans le fichier)
    vector<string> typeName;
    vector<vector<string>> typeAtoms;
    vector<int> typeCount;
    vector<int> headName;  // decalage de l'atome de tete (par type) d'apres la table de noms
    vector<int> headAuto;  // idem, detecte automatiquement
    vector<vector<float>> exposure;
    bool autoDone = false;
    vector<char> weakType;  // types a tete peu exposee (sterols, diacylglycerols, ceramides...)
    int nNoHead = 0;
};

static string trim(const char* s, int n) {
    int a = 0, b = n;
    while (a < b && (s[a] == ' ' || s[a] == '\0')) a++;
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\0')) b--;
    return string(s + a, b - a);
}

// resn / atn : tableaux de chaines de largeur fixe W
// molAtom (optionnel) : numero de molecule de chaque atome (p. ex. lu dans un fichier TPR) ;
//   sans lui, une molecule = suite d'atomes de meme numero et meme nom de residu.
// excl (optionnel) : noms de residus a exclure (une molecule en contenant un est ignoree).
static void build_topology(Topo& T, int natoms, const int* resid, const char* resn, const char* atn, int W,
                           const int* molAtom = nullptr, const std::set<string>* excl = nullptr) {
    std::set<string> nonlip;
    for (int k = 0; NONLIPID[k]; k++) nonlip.insert(NONLIPID[k]);
    T.natoms = natoms;
    std::map<string, int> typeId;
    int mol = 0, i = 0;
    while (i < natoms) {
        int j = i + 1;
        if (molAtom)
            while (j < natoms && molAtom[j] == molAtom[i]) j++;
        else
            while (j < natoms && resid[j] == resid[i] && strncmp(resn + (size_t)j * W, resn + (size_t)i * W, W) == 0) j++;
        mol++;
        string rn = trim(resn + (size_t)i * W, W);
        // une molecule est un lipide si au moins un de ses residus n'est ni solvant, ni ion, ni acide amine
        bool lipidRes = false, excluded = false;
        string prev;
        for (int a = i; a < j; a++) {
            const char* r = resn + (size_t)a * W;
            if (a > i && strncmp(r, resn + (size_t)(a - 1) * W, W) == 0) continue;
            string ra = trim(r, W);
            if (!nonlip.count(ra)) lipidRes = true;
            if (excl && excl->count(ra)) excluded = true;
        }
        if (lipidRes && !excluded && j - i >= 3) {
            string key = rn + "#" + std::to_string(j - i);
            auto it = typeId.find(key);
            int t;
            if (it == typeId.end()) {
                t = (int)T.typeName.size();
                typeId[key] = t;
                T.typeName.push_back(rn);
                vector<string> names;
                for (int a = i; a < j; a++) names.push_back(trim(atn + (size_t)a * W, W));
                int h = -1;
                for (int k = 0; HEADNAMES[k] && h < 0; k++)
                    for (int a = 0; a < j - i; a++)
                        if (names[a] == HEADNAMES[k]) {
                            h = a;
                            break;
                        }
                if (h < 0) {
                    h = 0;
                    T.nNoHead++;
                }
                T.headName.push_back(h);
                T.typeAtoms.push_back(names);
                T.typeCount.push_back(0);
            } else
                t = it->second;
            T.typeCount[t]++;
            T.start.push_back(i);
            T.len.push_back(j - i);
            T.type.push_back(t);
            T.molid.push_back(mol);
        }
        i = j;
    }
    T.nlip = (int)T.start.size();
}

// ----------------------------------------------------------------------------
// Recherche de voisins : grille de cellules en coordonnees fractionnaires,
// images periodiques explicites (valable en triclinique et pour petites boites)
// ----------------------------------------------------------------------------
struct Edge {
    int i, j;
    float dx, dy, dz, d2;  // r_j - r_i (image la plus proche), distance au carre
};

struct Grid {
    int k[3], ncell, n;
    vector<int> cellStart, order;
    vector<float> fs;  // coordonnees fractionnaires triees par cellule
    void build(const float* P, int n_, const Box& B, double rc) {
        n = n_;
        for (int d = 0; d < 3; d++) k[d] = std::max(1, std::min(512, (int)std::floor(B.w[d] / rc)));
        while ((long)k[0] * k[1] * k[2] > 8L * n + 64) {
            int m = 0;
            for (int d = 1; d < 3; d++)
                if (k[d] > k[m]) m = d;
            k[m] = std::max(1, k[m] / 2);
        }
        ncell = k[0] * k[1] * k[2];
        vector<int> cell(n);
        vector<float> f((size_t)3 * n);
        cellStart.assign(ncell + 1, 0);
        for (int i = 0; i < n; i++) {
            double q[3];
            B.frac(P + 3 * i, q);
            int c[3];
            for (int d = 0; d < 3; d++) {
                q[d] -= std::floor(q[d]);
                if (q[d] >= 1.0) q[d] = 0.0;
                c[d] = std::min(k[d] - 1, (int)(q[d] * k[d]));
                f[3 * i + d] = (float)q[d];
            }
            cell[i] = (c[0] * k[1] + c[1]) * k[2] + c[2];
            cellStart[cell[i] + 1]++;
        }
        for (int c = 0; c < ncell; c++) cellStart[c + 1] += cellStart[c];
        order.resize(n);
        fs.resize((size_t)3 * n);
        vector<int> pos(cellStart.begin(), cellStart.end() - 1);
        for (int i = 0; i < n; i++) {
            int p = pos[cell[i]]++;
            order[p] = i;
            fs[3 * p] = f[3 * i];
            fs[3 * p + 1] = f[3 * i + 1];
            fs[3 * p + 2] = f[3 * i + 2];
        }
    }
    // appelle fn(thread, a, b, d[3], d2) pour chaque paire (a<b) a distance < rc
    template <class F>
    void pairs(const Box& B, double rc, int nthreads, F&& fn) const {
        const float rc2 = (float)(rc * rc);
        const bool half = (k[0] >= 3 && k[1] >= 3 && k[2] >= 3);
        float h[3][3];
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) h[a][b] = (float)B.H[a][b];
#pragma omp parallel for schedule(dynamic, 8) num_threads(nthreads) if (nthreads > 1)
        for (int c = 0; c < ncell; c++) {
            if (cellStart[c] == cellStart[c + 1]) continue;
#ifdef _OPENMP
            int th = omp_get_thread_num();
#else
            int th = 0;
#endif
            int cz = c % k[2], cy = (c / k[2]) % k[1], cx = c / (k[1] * k[2]);
            if (half) {  // paires internes a la cellule
                for (int p = cellStart[c]; p < cellStart[c + 1]; p++)
                    for (int q = p + 1; q < cellStart[c + 1]; q++) {
                        const float ux = fs[3 * q] - fs[3 * p], uy = fs[3 * q + 1] - fs[3 * p + 1], uz = fs[3 * q + 2] - fs[3 * p + 2];
                        float d[3];
                        d[0] = ux * h[0][0] + uy * h[1][0] + uz * h[2][0];
                        d[1] = ux * h[0][1] + uy * h[1][1] + uz * h[2][1];
                        d[2] = ux * h[0][2] + uy * h[1][2] + uz * h[2][2];
                        const float d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        if (d2 < rc2) fn(th, order[p], order[q], d, d2);
                    }
            }
            for (int ox = -1; ox <= 1; ox++) {
                int nx = cx + ox;
                float sx = 0;
                if (nx < 0) { nx += k[0]; sx = -1; } else if (nx >= k[0]) { nx -= k[0]; sx = 1; }
                for (int oy = -1; oy <= 1; oy++) {
                    int ny = cy + oy;
                    float sy = 0;
                    if (ny < 0) { ny += k[1]; sy = -1; } else if (ny >= k[1]) { ny -= k[1]; sy = 1; }
                    for (int oz = -1; oz <= 1; oz++) {
                        int nz = cz + oz;
                        float sz = 0;
                        if (nz < 0) { nz += k[2]; sz = -1; } else if (nz >= k[2]) { nz -= k[2]; sz = 1; }
                        // demi-voisinage (13 cellules) quand la grille le permet : chaque paire vue une fois
                        if (half && (ox < 0 || (ox == 0 && (oy < 0 || (oy == 0 && oz <= 0))))) continue;
                        int c2 = (nx * k[1] + ny) * k[2] + nz;
                        for (int p = cellStart[c]; p < cellStart[c + 1]; p++) {
                            const int a = order[p];
                            const float fx = fs[3 * p] - sx, fy = fs[3 * p + 1] - sy, fz = fs[3 * p + 2] - sz;
                            for (int q = cellStart[c2]; q < cellStart[c2 + 1]; q++) {
                                const int b = order[q];
                                if (!half && b <= a) continue;
                                const float ux = fs[3 * q] - fx, uy = fs[3 * q + 1] - fy, uz = fs[3 * q + 2] - fz;
                                float d[3];
                                d[0] = ux * h[0][0] + uy * h[1][0] + uz * h[2][0];
                                d[1] = ux * h[0][1] + uy * h[1][1] + uz * h[2][1];
                                d[2] = ux * h[0][2] + uy * h[1][2] + uz * h[2][2];
                                const float d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                                if (d2 < rc2) fn(th, a, b, d, d2);
                            }
                        }
                    }
                }
            }
        }
    }
    // appelle fn(b, d2) pour tous les points b a distance < rc du point a (requete ponctuelle)
    template <class F>
    void around(const Box& B, const vector<int>& inv, int a, double rc, F&& fn) const {
        const float rc2 = (float)(rc * rc);
        const int p = inv[a];
        int c[3];
        for (int d = 0; d < 3; d++) c[d] = std::min(k[d] - 1, (int)(fs[3 * p + d] * k[d]));
        for (int ox = -1; ox <= 1; ox++)
            for (int oy = -1; oy <= 1; oy++)
                for (int oz = -1; oz <= 1; oz++) {
                    int o3[3] = {ox, oy, oz}, m[3];
                    float sh[3];
                    for (int d = 0; d < 3; d++) {
                        m[d] = c[d] + o3[d];
                        sh[d] = 0;
                        if (m[d] < 0) { m[d] += k[d]; sh[d] = -1; } else if (m[d] >= k[d]) { m[d] -= k[d]; sh[d] = 1; }
                    }
                    int c2 = (m[0] * k[1] + m[1]) * k[2] + m[2];
                    for (int q = cellStart[c2]; q < cellStart[c2 + 1]; q++) {
                        if (q == p) continue;
                        const double u[3] = {fs[3 * q] + sh[0] - fs[3 * p], fs[3 * q + 1] + sh[1] - fs[3 * p + 1],
                                             fs[3 * q + 2] + sh[2] - fs[3 * p + 2]};
                        float d2 = 0;
                        for (int d = 0; d < 3; d++) {
                            float x = (float)(u[0] * B.H[0][d] + u[1] * B.H[1][d] + u[2] * B.H[2][d]);
                            d2 += x * x;
                        }
                        if (d2 < rc2) fn(order[q], d2);
                    }
                }
    }
};

static void neighbor_edges(const float* P, int n, const Box& B, double rc, int nthreads, vector<Edge>& E) {
    Grid G;
    G.build(P, n, B, rc);
    int nt = std::max(1, nthreads);
    E.clear();
    if (nt == 1) {
        G.pairs(B, rc, 1, [&](int, int a, int b, const float* d, float d2) { E.push_back({a, b, d[0], d[1], d[2], d2}); });
        return;
    }
    vector<vector<Edge>> buf(nt);
    G.pairs(B, rc, nt, [&](int th, int a, int b, const float* d, float d2) { buf[th].push_back({a, b, d[0], d[1], d[2], d2}); });
    size_t tot = 0;
    for (auto& b : buf) tot += b.size();
    E.reserve(tot);
    for (auto& b : buf) E.insert(E.end(), b.begin(), b.end());
}

// ----------------------------------------------------------------------------
// Union-Find (compression de chemin + union par taille)
// ----------------------------------------------------------------------------
struct UF {
    vector<int> p, sz;
    void init(int n) {
        p.resize(n);
        std::iota(p.begin(), p.end(), 0);
        sz.assign(n, 1);
    }
    inline int find(int x) {
        while (p[x] != x) {
            p[x] = p[p[x]];
            x = p[x];
        }
        return x;
    }
    // renvoie la taille de la plus petite des deux composantes fusionnees (0 si deja unies)
    inline int unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return 0;
        if (sz[a] < sz[b]) std::swap(a, b);
        int s = sz[b];
        p[b] = a;
        sz[a] += s;
        return s;
    }
};

// ----------------------------------------------------------------------------
// Parametres et resultats
// ----------------------------------------------------------------------------
struct Params {
    int version = 4;
    double cutoff = 0;      // <=0 : automatique (v1,v2) ou valeur par defaut (v3,v4)
    double rmax = 2.5;      // v1/v2 : distance max exploree pour le seuil automatique (nm)
    double binw = 0.02;     // v1/v2 : resolution du balayage du seuil (nm)
    double sigfrac = 0.10;  // v1/v2 : une fusion est "significative" si la petite composante >= sigfrac*N
    double cdef = 1.5;      // v3/v4 : rayon de voisinage par defaut (nm)
    double rnorm = 2.0;     // v3 : rayon du voisinage pour l'ACP (nm)
    double costh = 0.5;     // v3/v4 : cosinus minimal entre normales voisines
    double qmin = 0.5;      // v4 : coherence minimale du voisinage pour un lipide de coeur
    double minfrac = 0.02;  // taille minimale d'un feuillet (fraction de N)
    int minsize = 10;       // taille minimale absolue d'un feuillet
    double rexp = 1.0;      // rayon pour le calcul d'exposition (detection des tetes)
    int attach = 1;         // v3/v4 : rattachement des orphelins : 0 = par orientation, 1 = plus proche voisin
    bool weakCore = true;   // v4 : les types a tete enfouie ne construisent pas les feuillets (rattaches ensuite)
    int nthreads = 1;
    bool namedHeads = false;  // v3/v4 : utiliser les tetes par nom plutot que detectees
};
struct Result {
    double cutoff = 0, t_total = 0, t_geom = 0, t_neigh = 0, t_graph = 0, t_heads = 0;
    int nleaflets = 0, nunassigned = 0, norphans = 0;
    long nedges = 0;
};

// ----------------------------------------------------------------------------
// Detection automatique de l'atome de reference de chaque type de lipide :
// atome le plus "expose", c.-a-d. le plus eloigne du barycentre des atomes
// lipidiques situes dans une sphere de rayon rexp autour de lui. Les atomes
// enfouis dans le coeur hydrophobe ont un environnement symetrique (expo ~ 0),
// ceux de l'interface ont un environnement hemispherique (expo maximale).
// Aucun nom d'atome ni solvant n'est necessaire. Calcule une seule fois.
// ----------------------------------------------------------------------------
static void detect_heads(Topo& T, const float* xyz, const Box& B, const Params& P) {
    int nt = (int)T.typeName.size();
    // tous les atomes lipidiques
    vector<int> lipOf, offOf;
    vector<float> X;
    for (int l = 0; l < T.nlip; l++)
        for (int a = 0; a < T.len[l]; a++) {
            lipOf.push_back(l);
            offOf.push_back(a);
            const float* r = xyz + 3 * (size_t)(T.start[l] + a);
            X.push_back(r[0]);
            X.push_back(r[1]);
            X.push_back(r[2]);
        }
    int na = (int)lipOf.size();
    // sous-echantillonnage : au plus maxMol molecules par type contribuent
    const int maxMol = 400;
    vector<char> use(T.nlip, 0);
    {
        vector<int> seen(nt, 0);
        for (int l = 0; l < T.nlip; l++) {
            int t = T.type[l];
            int stride = std::max(1, T.typeCount[t] / maxMol);
            if (seen[t] % stride == 0) use[l] = 1;
            seen[t]++;
        }
    }
    Grid G;
    G.build(X.data(), na, B, P.rexp);
    int nth = std::max(1, P.nthreads);
    vector<vector<float>> acc(nth, vector<float>((size_t)4 * na, 0.f));
    G.pairs(B, P.rexp, nth, [&](int th, int a, int b, const float* d, float) {
        if (lipOf[a] == lipOf[b]) return;  // seuls les atomes des autres molecules comptent
        float* v = acc[th].data();
        if (use[lipOf[a]]) { v[4 * a] += d[0]; v[4 * a + 1] += d[1]; v[4 * a + 2] += d[2]; v[4 * a + 3] += 1; }
        if (use[lipOf[b]]) { v[4 * b] -= d[0]; v[4 * b + 1] -= d[1]; v[4 * b + 2] -= d[2]; v[4 * b + 3] += 1; }
    });
    T.exposure.assign(nt, vector<float>());
    vector<vector<int>> cnt(nt);
    for (int t = 0; t < nt; t++) {
        T.exposure[t].assign(T.typeAtoms[t].size(), 0.f);
        cnt[t].assign(T.typeAtoms[t].size(), 0);
    }
    for (int a = 0; a < na; a++) {
        if (!use[lipOf[a]]) continue;
        float s[4] = {0, 0, 0, 0};
        for (int th = 0; th < nth; th++)
            for (int k = 0; k < 4; k++) s[k] += acc[th][4 * (size_t)a + k];
        float e = (float)P.rexp;  // atome isole : exposition maximale
        if (s[3] > 0) e = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]) / s[3];
        int t = T.type[lipOf[a]];
        T.exposure[t][offOf[a]] += e;
        cnt[t][offOf[a]]++;
    }
    T.headAuto.assign(nt, 0);
    for (int t = 0; t < nt; t++) {
        int best = 0;
        for (size_t a = 0; a < T.exposure[t].size(); a++) {
            if (cnt[t][a]) T.exposure[t][a] /= cnt[t][a];
            if (T.exposure[t][a] > T.exposure[t][best]) best = (int)a;
        }
        T.headAuto[t] = best;
    }
    // Seconde passe pour les types "enfouis" (exposition maximale faible : diacylglycerols,
    // ceramides, sterols...) dont le signal d'exposition est peu discriminant : on retient
    // l'atome en moyenne le plus proche des tetes des lipides "fiables" voisins.
    double wsum = 0, esum = 0;
    for (int t = 0; t < nt; t++) {
        wsum += T.typeCount[t];
        esum += (double)T.typeCount[t] * T.exposure[t][T.headAuto[t]];
    }
    const float ethr = (float)(0.6 * esum / std::max(1.0, wsum));
    vector<char> weak(nt, 0);
    int nweak = 0;
    for (int t = 0; t < nt; t++)
        if (T.exposure[t][T.headAuto[t]] < ethr) { weak[t] = 1; nweak++; }
    T.weakType = weak;
    if (nweak && nweak < nt) {
        // points : tetes fiables (indices < nh) puis atomes des lipides "faibles"
        vector<float> Y;
        vector<int> wl, wo;
        for (int l = 0; l < T.nlip; l++)
            if (!weak[T.type[l]]) {
                const float* r = xyz + 3 * (size_t)(T.start[l] + T.headAuto[T.type[l]]);
                Y.insert(Y.end(), r, r + 3);
            }
        int nh = (int)Y.size() / 3;
        for (int l = 0; l < T.nlip; l++)
            if (weak[T.type[l]] && use[l])
                for (int a = 0; a < T.len[l]; a++) {
                    const float* r = xyz + 3 * (size_t)(T.start[l] + a);
                    Y.insert(Y.end(), r, r + 3);
                    wl.push_back(l);
                    wo.push_back(a);
                }
        int ny = (int)Y.size() / 3;
        const double r2 = 2.0 * P.rexp;
        vector<float> dmin(ny - nh, (float)(r2 * r2));
        Grid G2;
        G2.build(Y.data(), ny, B, r2);
        G2.pairs(B, r2, 1, [&](int, int a, int b, const float*, float d2) {
            if (a < nh && b >= nh) { if (d2 < dmin[b - nh]) dmin[b - nh] = d2; }
        });
        vector<vector<double>> ds(nt);
        vector<vector<int>> dc(nt);
        for (int t = 0; t < nt; t++)
            if (weak[t]) { ds[t].assign(T.typeAtoms[t].size(), 0.0); dc[t].assign(T.typeAtoms[t].size(), 0); }
        for (int q = 0; q < ny - nh; q++) {
            int t = T.type[wl[q]];
            ds[t][wo[q]] += std::sqrt(dmin[q]);
            dc[t][wo[q]]++;
        }
        for (int t = 0; t < nt; t++) {
            if (!weak[t]) continue;
            int best = -1;
            for (size_t a = 0; a < ds[t].size(); a++)
                if (dc[t][a] && (best < 0 || ds[t][a] / dc[t][a] < ds[t][best] / dc[t][best])) best = (int)a;
            if (best >= 0) T.headAuto[t] = best;
        }
    }
    T.autoDone = true;
}

// ----------------------------------------------------------------------------
// Geometrie par lipide : position de l'atome de reference, barycentre et
// orientation (barycentre -> tete), molecule reconstruite a travers les PBC
// ----------------------------------------------------------------------------
static void lipid_geometry(const Topo& T, const float* xyz, const Box& B, const vector<int>& headOff, bool needOrient, int nthreads,
                           vector<float>& head, vector<float>& orient) {
    int n = T.nlip;
    head.resize((size_t)3 * n);
    if (needOrient) orient.resize((size_t)3 * n);
#pragma omp parallel for schedule(static) num_threads(nthreads) if (nthreads > 1)
    for (int l = 0; l < n; l++) {
        const float* r0 = xyz + 3 * (size_t)T.start[l];
        const int ho = headOff[T.type[l]];
        const float* rh = r0 + 3 * ho;
        head[3 * l] = rh[0];
        head[3 * l + 1] = rh[1];
        head[3 * l + 2] = rh[2];
        if (!needOrient) continue;
        // barycentre relatif a l'atome de tete, avec image minimale
        double s[3] = {0, 0, 0};
        for (int a = 0; a < T.len[l]; a++) {
            float d[3] = {r0[3 * a] - rh[0], r0[3 * a + 1] - rh[1], r0[3 * a + 2] - rh[2]};
            B.minimg(d);
            s[0] += d[0];
            s[1] += d[1];
            s[2] += d[2];
        }
        double nrm = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
        if (nrm < 1e-12) nrm = 1;
        orient[3 * l] = (float)(-s[0] / nrm);  // barycentre -> tete
        orient[3 * l + 1] = (float)(-s[1] / nrm);
        orient[3 * l + 2] = (float)(-s[2] / nrm);
    }
}

// ----------------------------------------------------------------------------
// Etiquetage final : composantes -> feuillets (tries par taille decroissante),
// petites composantes -> orphelins (-1) a rattacher
// ----------------------------------------------------------------------------
static int label_components(UF& uf, int n, const Params& P, vector<int>& lab) {
    int minsz = std::max(P.minsize, (int)std::ceil(P.minfrac * n));
    if (n < 2 * minsz) minsz = std::max(2, n / 4);
    vector<std::pair<int, int>> big;  // (taille, racine)
    for (int i = 0; i < n; i++)
        if (uf.p[i] == i && uf.sz[i] >= minsz) big.push_back({uf.sz[i], i});
    std::sort(big.begin(), big.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.first != b.first ? a.first > b.first : a.second < b.second;
    });
    vector<int> rootLab(n, -1);
    for (size_t k = 0; k < big.size(); k++) rootLab[big[k].second] = (int)k;
    lab.resize(n);
    for (int i = 0; i < n; i++) lab[i] = rootLab[uf.find(i)];
    return (int)big.size();
}

// rattachement des orphelins au feuillet du plus proche voisin deja assigne
// (propagation en quelques passes : un orphelin rattache peut servir de relais)
static void attach_nearest(const vector<Edge>& E, vector<int>& lab, int n) {
    for (int round = 0; round < 4; round++) {
        vector<float> best(n, 1e30f);
        vector<int> nl(n, -1);
        for (const Edge& e : E) {
            if (lab[e.i] < 0 && lab[e.j] >= 0 && e.d2 < best[e.i]) { best[e.i] = e.d2; nl[e.i] = lab[e.j]; }
            if (lab[e.j] < 0 && lab[e.i] >= 0 && e.d2 < best[e.j]) { best[e.j] = e.d2; nl[e.j] = lab[e.i]; }
        }
        int nnew = 0;
        for (int i = 0; i < n; i++)
            if (lab[i] < 0 && nl[i] >= 0) { lab[i] = nl[i]; nnew++; }
        if (!nnew) break;
    }
}

// dernier recours pour les rares orphelins sans voisin a distance c (p. ex. cholesterol au
// centre de la bicouche) : plus proche lipide assigne dans un rayon rfar (recherche directe,
// cout n_orphelins x N, negligeable). Au-dela, le lipide est hors membrane : non assigne.
static void attach_far(const vector<float>& head, const Box& B, double rfar, vector<int>& lab, int n) {
    vector<int> orph;
    for (int i = 0; i < n; i++)
        if (lab[i] < 0) orph.push_back(i);
    if (orph.empty()) return;
    Grid G;
    G.build(head.data(), n, B, rfar);
    vector<int> inv(n);
    for (int p = 0; p < n; p++) inv[G.order[p]] = p;
    vector<int> nl(orph.size(), -1);
    for (size_t k = 0; k < orph.size(); k++) {
        float best = 1e30f;
        G.around(B, inv, orph[k], rfar, [&](int j, float d2) {
            if (lab[j] >= 0 && d2 < best) { best = d2; nl[k] = lab[j]; }
        });
    }
    for (size_t k = 0; k < orph.size(); k++) lab[orph[k]] = nl[k];
}

// ----------------------------------------------------------------------------
// V1 / V2 : graphe + seuil automatique.
// On construit une seule fois les aretes jusqu'a rmax, on les range par classes
// de distance (tri par comptage) puis on les insere par distance croissante
// dans un Union-Find (equivalent a un clustering "single linkage"). La derniere
// fusion "significative" (deux grandes composantes) est la fusion des deux
// feuillets : le seuil est place au milieu de la plage ou les feuillets sont
// formes mais pas encore fusionnes.
// ----------------------------------------------------------------------------
static void graph_auto(const vector<float>& head, int n, const Box& B, const Params& P, vector<int>& lab, Result& R) {
    double t0 = now();
    static thread_local vector<Edge> E;  // tampon reutilise d'une trame a l'autre
    double rmax = P.cutoff > 0 ? P.cutoff : P.rmax;
    neighbor_edges(head.data(), n, B, rmax, P.nthreads, E);
    R.nedges = (long)E.size();
    double t1 = now();
    R.t_neigh += t1 - t0;
    double cut = P.cutoff;
    UF uf;
    if (cut <= 0) {
        int nb = (int)std::ceil(rmax / P.binw);
        vector<int> bstart(nb + 2, 0);
        vector<int> ebin(E.size());
        for (size_t e = 0; e < E.size(); e++) {
            int b = std::min(nb - 1, (int)(std::sqrt(E[e].d2) / P.binw));
            ebin[e] = b;
            bstart[b + 1]++;
        }
        for (int b = 0; b <= nb; b++) bstart[b + 1] += bstart[b];
        vector<int> sorted(E.size()), pos(bstart.begin(), bstart.end() - 1);
        for (size_t e = 0; e < E.size(); e++) sorted[pos[ebin[e]]++] = (int)e;
        uf.init(n);
        int sig = std::max(2, (int)std::ceil(P.sigfrac * n));
        vector<int> s1(nb, 0), s2(nb, 0);
        int lastSig = -1;
        bool changed = true;
        int c1 = 1, c2 = 1;
        for (int b = 0; b < nb; b++) {
            for (int q = bstart[b]; q < bstart[b + 1]; q++) {
                int s = uf.unite(E[sorted[q]].i, E[sorted[q]].j);
                if (s) changed = true;
                if (s >= sig) lastSig = b;
            }
            if (changed) {  // deux plus grandes composantes
                c1 = c2 = 0;
                for (int i = 0; i < n; i++)
                    if (uf.p[i] == i) {
                        int s = uf.sz[i];
                        if (s > c1) { c2 = c1; c1 = s; } else if (s > c2) c2 = s;
                    }
                changed = false;
            }
            s1[b] = c1;
            s2[b] = c2;
        }
        int hi;  // derniere classe ou les feuillets sont encore separes
        if (s2[nb - 1] >= sig || lastSig <= 0) hi = nb - 1;  // jamais fusionnes dans [0, rmax] : plateau
        else hi = lastSig - 1;
        int lo = hi;
        for (int b = 0; b <= hi; b++)
            if (s1[b] + s2[b] >= 0.9 * (s1[hi] + s2[hi]) && s2[b] >= 0.9 * s2[hi]) { lo = b; break; }
        cut = 0.5 * ((lo + 1) * P.binw + (hi + 1) * P.binw);
    }
    uf.init(n);
    float c2f = (float)(cut * cut);
    for (const Edge& e : E)
        if (e.d2 < c2f) uf.unite(e.i, e.j);
    R.cutoff = cut;
    R.nleaflets = label_components(uf, n, P, lab);
    for (int i = 0; i < n; i++) R.norphans += (lab[i] < 0);
    attach_nearest(E, lab, n);
    R.t_graph += now() - t1;
}

// plus petit vecteur propre d'une matrice symetrique 3x3 (Jacobi)
static void smallest_eigvec(double A[3][3], float* v) {
    double V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int sweep = 0; sweep < 12; sweep++) {
        double off = std::fabs(A[0][1]) + std::fabs(A[0][2]) + std::fabs(A[1][2]);
        if (off < 1e-14) break;
        for (int p = 0; p < 2; p++)
            for (int q = p + 1; q < 3; q++) {
                if (std::fabs(A[p][q]) < 1e-300) continue;
                double th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
                double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
                double c = 1 / std::sqrt(t * t + 1), s = t * c;
                for (int k = 0; k < 3; k++) {
                    double akp = A[k][p], akq = A[k][q];
                    A[k][p] = c * akp - s * akq;
                    A[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; k++) {
                    double apk = A[p][k], aqk = A[q][k];
                    A[p][k] = c * apk - s * aqk;
                    A[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; k++) {
                    double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = c * vkp - s * vkq;
                    V[k][q] = s * vkp + c * vkq;
                }
            }
    }
    int m = 0;
    for (int k = 1; k < 3; k++)
        if (A[k][k] < A[m][m]) m = k;
    v[0] = (float)V[0][m];
    v[1] = (float)V[1][m];
    v[2] = (float)V[2][m];
}

// ----------------------------------------------------------------------------
// V3 / V4 / V5 : graphe contraint par l'orientation
//   V3 : normale locale = ACP des tetes voisines de meme orientation (rayon
//        rnorm), signe donne par l'orientation moyenne ; arete si distance < c,
//        normales compatibles et deplacement quasi tangent.
//   V4 : normale locale = moyenne des orientations des voisins (rayon c) ;
//        arete si distance < c et normales compatibles.
//   V5 : pas de lissage : la "normale" est l'orientation brute du lipide.
// Les lipides restes isoles (orientation ambigue, p. ex. cholesterol au centre
// de la bicouche) sont ensuite rattaches au feuillet voisin dont la normale
// est la plus compatible avec leur propre orientation.
// ----------------------------------------------------------------------------
static void graph_oriented(const vector<float>& head, const vector<float>& orient, const vector<char>* lipWeak, int n,
                           const Box& B, const Params& P,
                           vector<int>& lab, Result& R) {
    double t0 = now();
    const int ver = P.version;
    double c = P.cutoff > 0 ? P.cutoff : P.cdef;
    double rs = (ver == 3) ? std::max(c, P.rnorm) : c;
    static thread_local vector<Edge> E;  // tampon reutilise d'une trame a l'autre
    neighbor_edges(head.data(), n, B, rs, P.nthreads, E);
    R.nedges = (long)E.size();
    double t1 = now();
    R.t_neigh += t1 - t0;
    const float c2 = (float)(c * c);
    const float* o = orient.data();
    vector<float> nv((size_t)3 * n);
    vector<char> core(n, 1);  // lipides "de coeur" : seuls ceux-ci construisent le graphe
    vector<int> cnt(n, 0);
    if (ver == 5) {
        nv = orient;
    } else if (ver == 4) {
        for (int i = 0; i < 3 * n; i++) nv[i] = o[i];
        for (const Edge& e : E) {
            if (e.d2 >= c2) continue;
            const int i = e.i, j = e.j;
            nv[3 * i] += o[3 * j]; nv[3 * i + 1] += o[3 * j + 1]; nv[3 * i + 2] += o[3 * j + 2];
            nv[3 * j] += o[3 * i]; nv[3 * j + 1] += o[3 * i + 1]; nv[3 * j + 2] += o[3 * i + 2];
        }
        for (const Edge& e : E)
            if (e.d2 < c2) { cnt[e.i]++; cnt[e.j]++; }
        for (int i = 0; i < n; i++) {
            float* v = &nv[3 * i];
            float nr = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            // coherence du voisinage q = |somme o_j| / N : proche de 1 dans un feuillet,
            // faible au centre de la bicouche (voisins des deux feuillets)
            float q = nr / (float)(cnt[i] + 1);
            if (nr < 1e-6f) { v[0] = o[3 * i]; v[1] = o[3 * i + 1]; v[2] = o[3 * i + 2]; } else { v[0] /= nr; v[1] /= nr; v[2] /= nr; }
            float self = o[3 * i] * v[0] + o[3 * i + 1] * v[1] + o[3 * i + 2] * v[2];
            core[i] = (q >= (float)P.qmin && self >= (float)P.costh) && !(lipWeak && (*lipWeak)[i]);
        }
    } else {  // ver == 3 : ACP
        const float rn2 = (float)(P.rnorm * P.rnorm);
        vector<double> S((size_t)10 * n, 0.0);  // somme d (3), somme d d^T (6), nombre
        vector<float> so(orient);               // somme des orientations
        for (const Edge& e : E) {
            if (e.d2 >= rn2) continue;
            const int i = e.i, j = e.j;
            if (o[3 * i] * o[3 * j] + o[3 * i + 1] * o[3 * j + 1] + o[3 * i + 2] * o[3 * j + 2] <= 0) continue;
            const double d[3] = {e.dx, e.dy, e.dz};
            double* a = &S[10 * (size_t)i];
            double* b = &S[10 * (size_t)j];
            a[0] += d[0]; a[1] += d[1]; a[2] += d[2];
            b[0] -= d[0]; b[1] -= d[1]; b[2] -= d[2];
            const double m[6] = {d[0] * d[0], d[0] * d[1], d[0] * d[2], d[1] * d[1], d[1] * d[2], d[2] * d[2]};
            for (int k = 0; k < 6; k++) { a[3 + k] += m[k]; b[3 + k] += m[k]; }
            a[9] += 1; b[9] += 1;
            for (int k = 0; k < 3; k++) { so[3 * i + k] += o[3 * j + k]; so[3 * j + k] += o[3 * i + k]; }
        }
#pragma omp parallel for schedule(static) num_threads(P.nthreads) if (P.nthreads > 1)
        for (int i = 0; i < n; i++) {
            const double* a = &S[10 * (size_t)i];
            float* v = &nv[3 * i];
            float* s = &so[3 * i];
            float sn = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
            if (a[9] < 3) {  // pas assez de voisins : orientation moyenne
                for (int k = 0; k < 3; k++) v[k] = sn > 0 ? s[k] / sn : o[3 * i + k];
                continue;
            }
            double m = a[9] + 1;  // le point lui-meme (d = 0) est inclus
            double mu[3] = {a[0] / m, a[1] / m, a[2] / m};
            double C[3][3];
            C[0][0] = a[3] / m - mu[0] * mu[0];
            C[0][1] = C[1][0] = a[4] / m - mu[0] * mu[1];
            C[0][2] = C[2][0] = a[5] / m - mu[0] * mu[2];
            C[1][1] = a[6] / m - mu[1] * mu[1];
            C[1][2] = C[2][1] = a[7] / m - mu[1] * mu[2];
            C[2][2] = a[8] / m - mu[2] * mu[2];
            smallest_eigvec(C, v);
            if (v[0] * s[0] + v[1] * s[1] + v[2] * s[2] < 0) { v[0] = -v[0]; v[1] = -v[1]; v[2] = -v[2]; }
        }
    }
    UF uf;
    uf.init(n);
    const float cth = (float)P.costh;
    for (const Edge& e : E) {
        if (e.d2 >= c2 || !core[e.i] || !core[e.j]) continue;
        const float* a = &nv[3 * e.i];
        const float* b = &nv[3 * e.j];
        if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] <= cth) continue;
        if (ver == 3) {  // deplacement quasi tangent a la surface (|d.n| < |d|/2 ou < 0.5 nm)
            float dn = 0.5f * std::fabs(e.dx * (a[0] + b[0]) + e.dy * (a[1] + b[1]) + e.dz * (a[2] + b[2]));
            if (dn > 0.5f && dn * dn > 0.25f * e.d2) continue;
        }
        uf.unite(e.i, e.j);
    }
    R.cutoff = c;
    R.nleaflets = label_components(uf, n, P, lab);
    // rattachement des orphelins : score = somme des (o_i . n_j) sur les voisins j de chaque feuillet
    int nl = R.nleaflets, norph = 0;
    vector<int> oid(n, -1);
    for (int i = 0; i < n; i++)
        if (lab[i] < 0) oid[i] = norph++;
    R.norphans = norph;
    if (P.attach == 1 && norph && nl) attach_nearest(E, lab, n);
    for (int round = 0; round < 4 && norph && nl && P.attach == 0; round++) {
        vector<float> score((size_t)norph * nl, 0.f);
        vector<char> has((size_t)norph * nl, 0);
        auto add = [&](int i, int j) {  // i orphelin, j assigne
            float s = o[3 * i] * nv[3 * j] + o[3 * i + 1] * nv[3 * j + 1] + o[3 * i + 2] * nv[3 * j + 2];
            score[(size_t)oid[i] * nl + lab[j]] += s;
            has[(size_t)oid[i] * nl + lab[j]] = 1;
        };
        for (const Edge& e : E) {
            if (e.d2 >= c2) continue;
            if (lab[e.i] < 0 && lab[e.j] >= 0) add(e.i, e.j);
            else if (lab[e.j] < 0 && lab[e.i] >= 0) add(e.j, e.i);
        }
        int nnew = 0;
        vector<int> newlab(norph, -1);
        for (int i = 0; i < n; i++) {
            if (oid[i] < 0 || lab[i] >= 0) continue;
            float best = -1e30f;
            for (int l = 0; l < nl; l++)
                if (has[(size_t)oid[i] * nl + l] && score[(size_t)oid[i] * nl + l] > best) {
                    best = score[(size_t)oid[i] * nl + l];
                    newlab[oid[i]] = l;
                }
        }
        for (int i = 0; i < n; i++)
            if (oid[i] >= 0 && lab[i] < 0 && newlab[oid[i]] >= 0) { lab[i] = newlab[oid[i]]; nnew++; }
        if (!nnew) break;
    }
    if (nl) attach_far(head, B, 2.0 * c, lab, n);
    R.t_graph += now() - t1;
}

// ----------------------------------------------------------------------------
// Analyse d'une trame
// ----------------------------------------------------------------------------
[[maybe_unused]] static void prepare_heads(Topo& T, const float* xyz, const Box& B, const Params& P, Result* R = nullptr) {
    if (T.autoDone) return;
    double t0 = now();
    detect_heads(T, xyz, B, P);
    if (R) R->t_heads = now() - t0;
}

static void analyse(const Topo& T, const float* xyz, const Box& B, const Params& P, vector<int>& lab, Result& R) {
    double t0 = now();
    bool named = (P.version == 1) || (P.namedHeads && P.version >= 3);
    const vector<int>& hoff = named ? T.headName : T.headAuto;
    vector<float> head, orient;
    lipid_geometry(T, xyz, B, hoff, P.version >= 3, P.nthreads, head, orient);
    R.t_geom = now() - t0;
    if (P.version <= 2) graph_auto(head, T.nlip, B, P, lab, R);
    else {
        // lipides dont le type a une tete enfouie : exclus de la construction des feuillets (V4)
        vector<char> lipWeak;
        bool useWeak = P.version == 4 && P.weakCore && !T.weakType.empty();
        if (useWeak) {
            lipWeak.resize(T.nlip);
            for (int i = 0; i < T.nlip; i++) lipWeak[i] = T.weakType[T.type[i]];
        }
        graph_oriented(head, orient, useWeak ? &lipWeak : nullptr, T.nlip, B, P, lab, R);
    }
    R.nunassigned = 0;
    for (int i = 0; i < T.nlip; i++) R.nunassigned += (lab[i] < 0);
    R.t_total = now() - t0;
}

// =============================================================================
// API C (bibliotheque partagee, utilisee depuis Python / ctypes)
// =============================================================================
extern "C" {
// molAtom : numero de molecule par atome (ou NULL) ; exclude : noms de residus separes par des virgules (ou NULL)
static std::set<string> split_csv(const char* s) {
    std::set<string> out;
    if (!s) return out;
    string cur;
    for (const char* p = s;; p++) {
        if (*p == ',' || *p == 0) {
            if (!cur.empty()) out.insert(cur);
            cur.clear();
            if (!*p) break;
        } else if (*p != ' ')
            cur += *p;
    }
    return out;
}
void* lf_new(int natoms, const int* resid, const char* resnames, const char* atomnames, int width, const int* molAtom,
             const char* exclude) {
    tune_malloc();
    Topo* T = new Topo;
    std::set<string> ex = split_csv(exclude);
    build_topology(*T, natoms, resid, resnames, atomnames, width, molAtom, ex.empty() ? nullptr : &ex);
    return T;
}
void lf_free(void* h) { delete (Topo*)h; }
int lf_nlip(void* h) { return ((Topo*)h)->nlip; }
void lf_molids(void* h, int* out) {
    Topo* T = (Topo*)h;
    for (int i = 0; i < T->nlip; i++) out[i] = T->molid[i];
}
// indice (0..natoms-1) de l'atome de reference de chaque lipide (autohead=1 : detecte)
void lf_head_atoms(void* h, int autohead, int* out) {
    Topo* T = (Topo*)h;
    for (int i = 0; i < T->nlip; i++) out[i] = T->start[i] + (autohead ? T->headAuto[T->type[i]] : T->headName[T->type[i]]);
}
static void box_from(const float* box9, Box& B) {
    double H[3][3];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) H[a][b] = box9[3 * a + b];
    B.set(H);
}
// detection des atomes de reference (a appeler une fois, sur la premiere trame)
double lf_prepare(void* h, const float* xyz, const float* box9, int nthreads) {
    Topo* T = (Topo*)h;
    Box B;
    box_from(box9, B);
    Params P;
    P.nthreads = nthreads;
    double t0 = now();
    T->autoDone = false;
    detect_heads(*T, xyz, B, P);
    return now() - t0;
}
// info[0]=seuil, [1]=temps total (s), [2]=nb feuillets, [3]=nb non assignes, [4]=nb orphelins rattaches ou non
int lf_run(void* h, const float* xyz, const float* box9, int version, double cutoff, int nthreads, int* labels, double* info) {
    Topo* T = (Topo*)h;
    Box B;
    box_from(box9, B);
    Params P;
    P.version = version;
    P.cutoff = cutoff;
    P.nthreads = nthreads;
    if (version >= 2 && !T->autoDone) return -1;
    vector<int> lab;
    Result R;
    analyse(*T, xyz, B, P, lab, R);
    for (int i = 0; i < T->nlip; i++) labels[i] = lab[i];
    if (info) {
        info[0] = R.cutoff;
        info[1] = R.t_total;
        info[2] = R.nleaflets;
        info[3] = R.nunassigned;
        info[4] = R.norphans;
    }
    return R.nleaflets;
}
}

// =============================================================================
// Programme principal : lecture GRO (mono ou multi-trames), analyse, sorties
// =============================================================================
#ifndef LF_LIB

struct GroFile {
    vector<char> buf;
    size_t pos = 0;
    int natoms = 0;
    vector<int> resid;
    vector<char> resn, atn;  // largeur 5
    vector<float> xyz;
    double H[3][3];
    string title;
    bool topoRead = false;
    bool load(const char* fn) {
        FILE* f = fopen(fn, "rb");
        if (!f) return false;
        fseek(f, 0, SEEK_END);
        long sz = ftell(f);
        fseek(f, 0, SEEK_SET);
        buf.resize(sz + 1);
        if (fread(buf.data(), 1, sz, f) != (size_t)sz) { fclose(f); return false; }
        buf[sz] = 0;
        fclose(f);
        return true;
    }
    inline const char* line(size_t& len) {
        if (pos >= buf.size() - 1) return nullptr;
        const char* s = buf.data() + pos;
        const char* e = (const char*)memchr(s, '\n', buf.size() - 1 - pos);
        if (!e) e = buf.data() + buf.size() - 1;
        len = e - s;
        pos = (e - buf.data()) + 1;
        if (len && s[len - 1] == '\r') len--;
        return s;
    }
    static inline float fastf(const char* s, int w) {
        int i = 0;
        while (i < w && s[i] == ' ') i++;
        bool neg = false;
        if (i < w && s[i] == '-') { neg = true; i++; }
        double v = 0;
        while (i < w && s[i] >= '0' && s[i] <= '9') v = v * 10 + (s[i++] - '0');
        if (i < w && s[i] == '.') {
            i++;
            double sc = 0.1;
            while (i < w && s[i] >= '0' && s[i] <= '9') { v += (s[i++] - '0') * sc; sc *= 0.1; }
        }
        return (float)(neg ? -v : v);
    }
    // lit la trame suivante ; renvoie false en fin de fichier
    bool next() {
        size_t len;
        const char* s = line(len);
        if (!s) return false;
        title.assign(s, len);
        s = line(len);
        if (!s) return false;
        int n = atoi(s);
        if (n <= 0) return false;
        if (!topoRead) {
            natoms = n;
            resid.resize(n);
            resn.assign((size_t)5 * n, ' ');
            atn.assign((size_t)5 * n, ' ');
            xyz.resize((size_t)3 * n);
        } else if (n != natoms) {
            fprintf(stderr, "Erreur : nombre d'atomes variable entre trames\n");
            return false;
        }
        int w = 8;
        for (int i = 0; i < n; i++) {
            s = line(len);
            if (!s || len < 39) { fprintf(stderr, "Erreur : ligne %d trop courte\n", i + 3); return false; }
            if (i == 0) {  // largeur des champs de coordonnees = distance entre deux points decimaux
                const char* p1 = (const char*)memchr(s + 20, '.', len - 20);
                const char* p2 = p1 ? (const char*)memchr(p1 + 1, '.', len - (p1 + 1 - s)) : nullptr;
                if (p1 && p2) w = (int)(p2 - p1);
            }
            if (!topoRead) {
                char tmp[6];
                memcpy(tmp, s, 5);
                tmp[5] = 0;
                resid[i] = atoi(tmp);
                memcpy(&resn[(size_t)5 * i], s + 5, 5);
                memcpy(&atn[(size_t)5 * i], s + 10, 5);
            }
            xyz[3 * (size_t)i] = fastf(s + 20, w);
            xyz[3 * (size_t)i + 1] = fastf(s + 20 + w, w);
            xyz[3 * (size_t)i + 2] = fastf(s + 20 + 2 * w, w);
        }
        s = line(len);
        if (!s) return false;
        double b[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        string bl(s, len);
        sscanf(bl.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf", b, b + 1, b + 2, b + 3, b + 4, b + 5, b + 6, b + 7, b + 8);
        H[0][0] = b[0]; H[1][1] = b[1]; H[2][2] = b[2];
        H[0][1] = b[3]; H[0][2] = b[4]; H[1][0] = b[5];
        H[1][2] = b[6]; H[2][0] = b[7]; H[2][1] = b[8];
        topoRead = true;
        return true;
    }
};

static vector<vector<int>> read_ref(const char* fn) {
    vector<vector<int>> g;
    FILE* f = fopen(fn, "r");
    if (!f) return g;
    char l[4096];
    while (fgets(l, sizeof(l), f)) {
        if (strchr(l, '[')) { g.push_back({}); continue; }
        if (g.empty()) continue;
        char* p = l;
        char* e;
        for (;;) {
            long v = strtol(p, &e, 10);
            if (e == p) break;
            g.back().push_back((int)v);
            p = e;
        }
    }
    fclose(f);
    return g;
}

// accord avec la reference : meilleure correspondance gloutonne feuillet trouve <-> feuillet de reference
static double accuracy(const Topo& T, const vector<int>& lab, int nl, const vector<vector<int>>& ref, int& nref, int& nmiss) {
    std::map<int, int> refOf;
    for (size_t k = 0; k < ref.size(); k++)
        for (int id : ref[k]) refOf[id] = (int)k;
    int nr = (int)ref.size();
    vector<vector<int>> M(std::max(nl, 1), vector<int>(nr, 0));
    nref = 0;
    nmiss = 0;
    for (int i = 0; i < T.nlip; i++) {
        auto it = refOf.find(T.molid[i]);
        if (it == refOf.end()) continue;
        nref++;
        if (lab[i] < 0) { nmiss++; continue; }
        M[lab[i]][it->second]++;
    }
    vector<char> used(nr, 0);
    int ok = 0;
    for (int l = 0; l < nl; l++) {
        int best = -1;
        for (int k = 0; k < nr; k++)
            if (!used[k] && (best < 0 || M[l][k] > M[l][best])) best = k;
        if (best < 0) break;
        used[best] = 1;
        ok += M[l][best];
    }
    return nref ? (double)ok / nref : 0;
}

static void write_leaflets(const char* fn, const Topo& T, const vector<int>& lab, int nl) {
    FILE* f = fopen(fn, "w");
    if (!f) { fprintf(stderr, "Impossible d'ecrire %s\n", fn); return; }
    for (int l = -1; l < nl; l++) {
        vector<int> ids;
        for (int i = 0; i < T.nlip; i++)
            if (lab[i] == l) ids.push_back(T.molid[i]);
        if (l < 0 && ids.empty()) continue;
        if (l < 0) fprintf(f, "[ non_assignes ]\n");
        else fprintf(f, "[ membrane_1_leaflet_%d ]\n", l + 1);
        for (size_t k = 0; k < ids.size(); k++) fprintf(f, "%6d%s", ids[k], (k % 15 == 14 || k + 1 == ids.size()) ? "\n" : "");
    }
    fclose(f);
}

static void usage() {
    printf("Usage : leaflet -f structure.gro [options]\n"
           "  -m N        version de l'algorithme : 1, 2, 3, 4 (defaut), 5\n"
           "  -o fichier  fichier de sortie des feuillets (format leaflets.txt)\n"
           "  -r fichier  fichier de reference (calcul du taux d'accord)\n"
           "  -c seuil    seuil de distance en nm (defaut : automatique en v1/v2, 1.5 en v3/v4)\n"
           "  -t N        nombre de threads OpenMP (defaut 1)\n"
           "  -n N        repeter l'analyse N fois (mesure du temps, defaut 1)\n"
           "  -P          trajectoire : analyser les trames en parallele (N threads sur les trames)\n"
           "  -s fichier  trajectoire : resume par trame au format CSV\n"
           "  --named     v3/v4 : atomes de reference par nom plutot que detectes\n"
           "  --molsizes fichier  decoupage en molecules (lignes \"nb_molecules atomes_par_molecule\", cf. scripts/tpr_molsizes.py)\n"
           "  -x RES1,RES2  exclure les molecules contenant ces residus (p. ex. -x S14 pour la surfactine)\n"
           "  --heads     afficher l'atome de reference detecte pour chaque type de lipide\n"
           "  -q          sortie minimale (une ligne de resume lisible par script)\n");
}

int main(int argc, char** argv) {
    const char *fn = nullptr, *out = nullptr, *reffn = nullptr, *csv = nullptr, *molfn = nullptr, *exclude = nullptr;
    Params P;
    int nrep = 1;
    bool quiet = false, showHeads = false, parFrames = false;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "-f" && i + 1 < argc) fn = argv[++i];
        else if (a == "-o" && i + 1 < argc) out = argv[++i];
        else if (a == "-r" && i + 1 < argc) reffn = argv[++i];
        else if (a == "-s" && i + 1 < argc) csv = argv[++i];
        else if (a == "-m" && i + 1 < argc) P.version = atoi(argv[++i]);
        else if (a == "-c" && i + 1 < argc) P.cutoff = atof(argv[++i]);
        else if (a == "-t" && i + 1 < argc) P.nthreads = atoi(argv[++i]);
        else if (a == "-n" && i + 1 < argc) nrep = atoi(argv[++i]);
        else if (a == "--costh" && i + 1 < argc) P.costh = atof(argv[++i]);
        else if (a == "--qmin" && i + 1 < argc) P.qmin = atof(argv[++i]);
        else if (a == "--rmax" && i + 1 < argc) P.rmax = atof(argv[++i]);
        else if (a == "--rexp" && i + 1 < argc) P.rexp = atof(argv[++i]);
        else if (a == "--attach" && i + 1 < argc) P.attach = atoi(argv[++i]);
        else if (a == "--molsizes" && i + 1 < argc) molfn = argv[++i];
        else if (a == "-x" && i + 1 < argc) exclude = argv[++i];
        else if (a == "--noweak") P.weakCore = false;
        else if (a == "--named") P.namedHeads = true;
        else if (a == "--heads") showHeads = true;
        else if (a == "-P") parFrames = true;
        else if (a == "-q") quiet = true;
        else { usage(); return a == "-h" ? 0 : 1; }
    }
    if (!fn || P.version < 1 || P.version > 5) { usage(); return 1; }
    tune_malloc();
    double t0 = now();
    GroFile G;
    if (!G.load(fn) || !G.next()) { fprintf(stderr, "Erreur de lecture de %s\n", fn); return 1; }
    double tread = now() - t0;
    t0 = now();
    Topo T;
    // decoupage en molecules : par residu (defaut) ou d'apres un fichier "nombre_de_molecules atomes_par_molecule"
    vector<int> molAtom;
    if (molfn) {
        FILE* fm = fopen(molfn, "r");
        if (!fm) { fprintf(stderr, "Impossible de lire %s\n", molfn); return 1; }
        long nm, na;
        int id = 0;
        while (fscanf(fm, "%ld %ld", &nm, &na) == 2)
            for (long m = 0; m < nm; m++, id++)
                for (long a = 0; a < na; a++) molAtom.push_back(id);
        fclose(fm);
        if ((int)molAtom.size() != G.natoms) {
            fprintf(stderr, "Erreur : %s decrit %zu atomes, le fichier GRO en contient %d\n", molfn, molAtom.size(), G.natoms);
            return 1;
        }
    }
    std::set<string> ex = split_csv(exclude);
    build_topology(T, G.natoms, G.resid.data(), G.resn.data(), G.atn.data(), 5, molfn ? molAtom.data() : nullptr,
                   ex.empty() ? nullptr : &ex);
    double ttopo = now() - t0;
    if (T.nlip == 0) { fprintf(stderr, "Aucun lipide trouve\n"); return 1; }
    Box B;
    B.set(G.H);
    Result R0;
    bool needAuto = (P.version == 2) || (P.version >= 3 && !P.namedHeads);
    if (needAuto || showHeads) prepare_heads(T, G.xyz.data(), B, P, &R0);
    if (showHeads) {
        printf("# type  nb_molecules  tete(nom)  tete(auto)  exposition_max\n");
        for (size_t t = 0; t < T.typeName.size(); t++)
            printf("%-6s %6d  %-5s  %-5s  %.3f\n", T.typeName[t].c_str(), T.typeCount[t], T.typeAtoms[t][T.headName[t]].c_str(),
                   T.typeAtoms[t][T.headAuto[t]].c_str(), T.exposure[t][T.headAuto[t]]);
    }
    vector<vector<int>> ref;
    if (reffn) ref = read_ref(reffn);

    // ---- premiere trame (eventuellement repetee pour la mesure du temps)
    vector<int> lab;
    Result R;
    double tbest = 1e30, tsum = 0;
    for (int r = 0; r < nrep; r++) {
        R = Result();
        analyse(T, G.xyz.data(), B, P, lab, R);
        tbest = std::min(tbest, R.t_total);
        tsum += R.t_total;
    }
    vector<int> sizes(R.nleaflets, 0);
    for (int i = 0; i < T.nlip; i++)
        if (lab[i] >= 0) sizes[lab[i]]++;
    double acc = -1;
    int nref = 0, nmiss = 0;
    if (!ref.empty()) acc = accuracy(T, lab, R.nleaflets, ref, nref, nmiss);
    if (quiet) {
        printf("version=%d nlip=%d nleaflets=%d unassigned=%d orphans=%d cutoff=%.3f edges=%ld t_best_ms=%.4f t_mean_ms=%.4f "
               "t_heads_ms=%.3f t_read_ms=%.2f acc=%.6f nref=%d sizes=",
               P.version, T.nlip, R.nleaflets, R.nunassigned, R.norphans, R.cutoff, R.nedges, tbest * 1e3, tsum / nrep * 1e3,
               R0.t_heads * 1e3, tread * 1e3, acc, nref);
        for (size_t k = 0; k < sizes.size(); k++) printf("%d%s", sizes[k], k + 1 < sizes.size() ? "," : "");
        printf("\n");
    } else {
        printf("Fichier            : %s\n", fn);
        printf("Atomes / lipides   : %d / %d (%zu types)\n", G.natoms, T.nlip, T.typeName.size());
        printf("Version            : %d   (threads : %d)\n", P.version, P.nthreads);
        printf("Seuil utilise      : %.3f nm   (%ld aretes candidates)\n", R.cutoff, R.nedges);
        printf("Feuillets trouves  : %d  -> tailles :", R.nleaflets);
        for (int s : sizes) printf(" %d", s);
        printf("\nOrphelins rattaches: %d   non assignes : %d\n", R.norphans - R.nunassigned, R.nunassigned);
        if (acc >= 0) printf("Accord reference   : %.4f %%  (%d lipides de reference, %d non assignes)\n", 100 * acc, nref, nmiss);
        printf("Temps lecture GRO  : %.3f ms   topologie : %.3f ms\n", tread * 1e3, ttopo * 1e3);
        if (needAuto) printf("Detection des tetes: %.3f ms (une seule fois)\n", R0.t_heads * 1e3);
        printf("Temps analyse      : %.4f ms (meilleur de %d ; geometrie %.4f, voisins %.4f, graphe %.4f)\n", tbest * 1e3, nrep,
               R.t_geom * 1e3, R.t_neigh * 1e3, R.t_graph * 1e3);
    }
    if (out) write_leaflets(out, T, lab, R.nleaflets);

    // ---- trames suivantes (fichier GRO multi-trames = trajectoire)
    vector<vector<float>> frames;
    vector<Box> boxes;
    double tr0 = now();
    while (G.next()) {
        frames.push_back(G.xyz);
        Box b;
        b.set(G.H);
        boxes.push_back(b);
    }
    double treadtraj = now() - tr0;
    if (!frames.empty()) {
        int nf = (int)frames.size();
        vector<vector<int>> labs(nf);
        vector<Result> res(nf);
        Params Pf = P;
        int nth = P.nthreads;
        if (parFrames) Pf.nthreads = 1;  // parallelisme sur les trames plutot que dans la trame
        double ta = now();
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth) if (parFrames && nth > 1)
        for (int f = 0; f < nf; f++) analyse(T, frames[f].data(), boxes[f], Pf, labs[f], res[f]);
        double tanal = now() - ta;
        // coherence des etiquettes d'une trame a l'autre + comptage des changements de feuillet
        FILE* fc = csv ? fopen(csv, "w") : nullptr;
        if (fc) fprintf(fc, "trame,n_feuillets,taille_1,taille_2,non_assignes,changements_vs_trame0,temps_ms\n");
        long totchg = 0;
        for (int f = 0; f < nf; f++) {
            int nl = res[f].nleaflets;
            // matrice de recouvrement avec la trame 0
            vector<vector<int>> M(std::max(nl, 1), vector<int>(std::max(R.nleaflets, 1), 0));
            for (int i = 0; i < T.nlip; i++)
                if (labs[f][i] >= 0 && lab[i] >= 0) M[labs[f][i]][lab[i]]++;
            vector<int> map(std::max(nl, 1), -1);
            vector<char> used(std::max(R.nleaflets, 1), 0);
            for (int l = 0; l < nl; l++) {
                int best = -1;
                for (int k = 0; k < R.nleaflets; k++)
                    if (!used[k] && (best < 0 || M[l][k] > M[l][best])) best = k;
                if (best >= 0) { used[best] = 1; map[l] = best; }
            }
            int chg = 0, s[2] = {0, 0};
            for (int i = 0; i < T.nlip; i++) {
                int l = labs[f][i] >= 0 ? map[labs[f][i]] : -1;
                labs[f][i] = l;
                if (l != lab[i]) chg++;
                if (l == 0 || l == 1) s[l]++;
            }
            totchg += chg;
            if (fc) fprintf(fc, "%d,%d,%d,%d,%d,%d,%.4f\n", f + 1, nl, s[0], s[1], res[f].nunassigned, chg, res[f].t_total * 1e3);
        }
        if (fc) fclose(fc);
        printf("Trajectoire        : %d trames supplementaires, lecture %.1f ms, analyse %.2f ms (%.4f ms/trame, %.0f trames/s)%s\n", nf,
               treadtraj * 1e3, tanal * 1e3, tanal * 1e3 / nf, nf / tanal, parFrames ? " [trames en parallele]" : "");
        printf("Changements de feuillet cumules par rapport a la trame 0 : %ld\n", totchg);
    }
    return 0;
}
#endif
