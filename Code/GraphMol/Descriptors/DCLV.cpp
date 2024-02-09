#define _USE_MATH_DEFINES

#include <iostream>
#include <limits>
#include <string>
#include <list>
#include <cmath>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include "DCLV.h"

#define VOXORDER 16
#define CHECKMAX 20
#define ATOMPOOL 16
#define HetAtmFlag 0x01
#define WaterFlag 0x10

namespace RDKit {
namespace Descriptors {

struct CheckType {
  const char* name;
  const unsigned int count;
};

static const CheckType residueCheck[CHECKMAX] = {
    {"ALA", 5},  /*  0 */
    {"GLY", 4},  /*  1 */
    {"LEU", 8},  /*  2 */
    {"SER", 6},  /*  3 */
    {"VAL", 7},  /*  4 */
    {"THR", 7},  /*  5 */
    {"LYS", 9},  /*  6 */
    {"ASP", 8},  /*  7 */
    {"ILE", 8},  /*  8 */
    {"ASN", 8},  /*  9 */
    {"GLU", 9},  /* 10 */
    {"PRO", 7},  /* 11 */
    {"ARG", 11}, /* 12 */
    {"PHE", 11}, /* 13 */
    {"GLN", 9},  /* 14 */
    {"TYR", 12}, /* 15 */
    {"HIS", 10}, /* 16 */
    {"CYS", 6},  /* 17 */
    {"MET", 8},  /* 18 */
    {"TRP", 14}  /* 19 */
};

struct DotStruct {
  double v[3];
  double area;
};

struct ElemStruct {
  std::vector<DotStruct> dots;
  double radius2;
  double radius;
  long count;
};

class AtomRecord {
 public:
  float radius; /* VdW */
  unsigned short atmSerNo;
  unsigned short resSerNo;
  float bFactor;
  float x, y, z;
  std::string atmName;
  std::string resName;
  std::string insert;
  std::string chain;
  int hetAtmFlag;
  bool solventFlag;
  char flag;
  ElemStruct* elem;

  bool isSolvent() {
    switch (resName[0]) {
      case 'D':
        return resName == "DOD" || resName == "D20";

      case 'H':
        return resName == "HOH" || resName == "H20";

      case 'S':
        return resName == "SOL" || resName == "SO4" || resName == "SUL";

      case 'W':
        return resName == "WAT";
      case 'T':
        return resName == "TIP";
      case 'P':
        return resName == "P04";
      default:
        return false;
    }
  }

  void initFlag() {
    flag = 0;  // initialise flag
    if (isSolvent()) {
      flag |= HetAtmFlag | WaterFlag;
    } else if (hetAtmFlag) {
      flag |= HetAtmFlag;
    }
  };

  // constructor to populate record
  AtomRecord(const Atom* atm, const Conformer cnf) {
    const AtomMonomerInfo* info = atm->getMonomerInfo();
    unsigned int i;
    solventFlag = false;

    if (info) {
      atmName = info->getName();
      resName = ((AtomPDBResidueInfo*)info)->getResidueName();

      if (isSolvent()) {
        solventFlag = true;
        return;
      }

      atmSerNo = ((AtomPDBResidueInfo*)info)->getSerialNumber();
      resSerNo = ((AtomPDBResidueInfo*)info)->getResidueNumber();
      x = cnf.getAtomPos(atm->getIdx()).x;
      y = cnf.getAtomPos(atm->getIdx()).y;
      z = cnf.getAtomPos(atm->getIdx()).z;
      insert = ((AtomPDBResidueInfo*)info)->getInsertionCode();
      chain = ((AtomPDBResidueInfo*)info)->getChainId();

      hetAtmFlag = 1;
      for (i = 0; i < CHECKMAX; i++) {
        if (residueCheck[i].name == resName) {
          hetAtmFlag = 0;
        }
      }
    } else {  // Molecule not from PDB (Ligand only)
      atmName = " " + atm->getSymbol();
      resName = "UNK";

      if (isSolvent()) {
        solventFlag = true;
        return;
      }

      atmSerNo = 1 + atm->getIdx();
      resSerNo = 0;
      x = cnf.getAtomPos(atm->getIdx()).x;
      y = cnf.getAtomPos(atm->getIdx()).y;
      z = cnf.getAtomPos(atm->getIdx()).z;
      hetAtmFlag = 1;
    }
  }
};

struct AtomList {
  const AtomRecord* ptr[ATOMPOOL];
  AtomList* next;
  unsigned int count;

  AtomList() : next(nullptr) {}
};

static void normalise(double* v) {
  double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

static double triangleArea(const double* p, const double* q, const double* r) {
  double a[3], b[3], c[3];

  a[0] = q[0] - p[0];
  b[0] = r[0] - p[0];
  a[1] = q[1] - p[1];
  b[1] = r[1] - p[1];
  a[2] = q[2] - p[2];
  b[2] = r[2] - p[2];

  c[0] = a[1] * b[2] - b[1] * a[2];
  c[1] = a[2] * b[0] - b[2] * a[0];
  c[2] = a[0] * b[1] - b[0] * a[1];

  double len = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
  return (0.5 * len);
}

static int within(const AtomRecord& src, const AtomRecord* dst, double dist) {
  double dx = src.x - dst->x;
  double dy = src.y - dst->y;
  double dz = src.z - dst->z;
  return ((dx * dx + dy * dy + dz * dz) < (dist * dist));
}

static void checkResidue(const AtomRecord* ptr, unsigned int count) {
  for (unsigned int i = 0; i < CHECKMAX; i++) {
    if (ptr->resName.compare(0, 3, residueCheck[i].name, 0, 3) == 0) {
      if (residueCheck[i].count != count) {
        std::cout << "Warning: Atom Count for residue " + ptr->resName +
                         std::to_string(ptr->resSerNo) + ptr->chain +
                         " does not match expected count"
                  << std::endl;
      }
      return;
    }
  }
}

static bool sameResidue(const AtomRecord* ptr1, const AtomRecord* ptr2) {
  if ((ptr1->flag | ptr2->flag) & HetAtmFlag) return (false);

  return ((ptr1->resSerNo == ptr2->resSerNo) &&
          (ptr1->resName.compare(0, 3, ptr2->resName, 0, 3) == 0) &&
          (ptr1->insert == ptr2->insert) && (ptr1->chain == ptr2->chain));
}

static double calculateCompactness(double d_surfaceArea, double d_totalVolume) {
  return d_surfaceArea / cbrt(36.0 * M_PI * d_totalVolume * d_totalVolume);
};

struct State {
  ElemStruct elemA;  // Alpha Carbon               //
  ElemStruct elemC;  // Peptide Carbon    Carbon   //
  ElemStruct elemN;  // Peptide Nitrogen  Nitrogen //
  ElemStruct elemO;  // Peptide Oxygen    Oxygen   //
  ElemStruct elemS;  //                   Sulphur  //
  ElemStruct elemX;  // Sidechain Atom    Unknown  //

  AtomList* grid[VOXORDER][VOXORDER][VOXORDER];
  AtomList* freeList;

  AtomList* neighbours;
  const AtomRecord* recordCache;

  unsigned long totalDots;
  ElemStruct standardDots;
  double standardArea;

  double voxX;
  double voxY;
  double voxZ;
  double voxU;
  double voxV;
  double voxW;
  double cenX;
  double cenY;
  double cenZ;

  void tesselate(const double* p, const double* q, const double* r,
                 unsigned int d) {
    double u[3], v[3], w[3];
    double* n;

    if (d--) {
      u[0] = p[0] + q[0];
      v[0] = q[0] + r[0];
      w[0] = r[0] + p[0];
      u[1] = p[1] + q[1];
      v[1] = q[1] + r[1];
      w[1] = r[1] + p[1];
      u[2] = p[2] + q[2];
      v[2] = q[2] + r[2];
      w[2] = r[2] + p[2];
      normalise(u);
      normalise(v);
      normalise(w);

      tesselate(u, v, w, d);
      tesselate(p, u, w, d);
      tesselate(u, q, v, d);
      tesselate(w, v, r, d);
    } else {
      double area = triangleArea(p, q, r);
      standardDots.dots[standardDots.count].area = area;
      n = standardDots.dots[standardDots.count++].v;
      standardArea += area;

      n[0] = p[0] + q[0] + r[0];
      n[1] = p[1] + q[1] + r[1];
      n[2] = p[2] + q[2] + r[2];
      normalise(n);
    }
  }

  void generateElemPoints(ElemStruct* elem, double rad, double probeRadius,
                          int dotDensity) {
    double x, y, z, p, q, xy;
    unsigned int vert;

    elem->radius = rad;
    rad = rad + probeRadius;
    elem->radius2 = rad * rad;

    if (dotDensity) {
      long count = ((4.0 * M_PI) * rad * rad * dotDensity);
      std::vector<DotStruct> dots(count);

      unsigned int equat = sqrt(M_PI * count);
      if (!(vert = equat >> 1)) vert = 1;

      unsigned int i = 0;
      for (unsigned int j = 0; (i < count) && (j < vert); j++) {
        p = (M_PI * j) / (double)vert;
        z = cos(p);
        xy = sin(p);
        unsigned int horz = (equat * xy);
        if (!horz) {
          horz = 1;
        }

        for (unsigned int k = 0; (i < count) && (k < horz); k++) {
          q = (2.0 * M_PI * k) / (double)horz;
          x = xy * sin(q);
          y = xy * cos(q);

          dots[i].v[0] = rad * x;
          dots[i].v[1] = rad * y;
          dots[i].v[2] = rad * z;
          i++;
        }
      }

      count = i;
      double area = ((4.0 * M_PI) * elem->radius2) / count;
      for (DotStruct& dot : dots) {
        dot.area = area;
      }

      elem->count = count;
    } else {
      elem->count = standardDots.count;
      elem->dots = standardDots.dots;
    }
  }

  void insertAtomList(AtomList** list, const AtomRecord* ptr) {
    while (*list && ((*list)->count == ATOMPOOL)) {
      list = &(*list)->next;
    }

    if (!(*list)) {
      if (freeList) {
        *list = freeList;
        freeList = (*list)->next;
        (*list)->next = nullptr;
        (*list)->count = 0;
      } else {
        (*list) = (AtomList*)malloc(sizeof(AtomList));
        (*list)->next = nullptr;
        (*list)->count = 0;
      }
    }
    (*list)->ptr[(*list)->count++] = ptr;
  }

  void freeAtomList(AtomList* ptr) {
    while (ptr) {
      AtomList* next = ptr->next;
      ptr->next = freeList;
      freeList = ptr;
      ptr = next;
    }
  }

  ElemStruct* getAtomElem(std::string atmName, bool typeFlag) {
    const char* name = atmName.c_str();

    if (typeFlag) { /* Only Recognise Backbone Atoms! */
      if ((name[1] == 'C') && (name[2] == 'A')) {
        return &elemA;
      } else if (name[2] == ' ')
        switch (name[1]) {
          case 'C':
            return &elemC;
          case 'N':
            return &elemN;
          case 'O':
            return &elemO;
        }
    } else { /* RasMol */
      switch (name[1]) {
        case 'C':
          return &elemC;
        case 'N':
          return &elemN;
        case 'O':
          return &elemO;
        case 'S':
          return &elemS;
      }
    }
    return (&elemX);
  }

  AtomList* findNeighbours(const AtomRecord& atom, double range) {
    double maxRadius = 1.87;
    double maxDist = range + maxRadius;

    int lx = voxX * (atom.x - maxDist - voxU);
    if (lx < 0) {
      lx = 0;
    }
    int ly = voxY * (atom.y - maxDist - voxV);
    if (ly < 0) {
      ly = 0;
    }
    int lz = voxZ * (atom.z - maxDist - voxW);
    if (lz < 0) {
      lz = 0;
    }

    int ux = voxX * (atom.x + maxDist - voxU);
    if (ux >= VOXORDER) {
      ux = VOXORDER - 1;
    }
    int uy = voxY * (atom.y + maxDist - voxV);
    if (uy >= VOXORDER) {
      uy = VOXORDER - 1;
    }
    int uz = voxZ * (atom.z + maxDist - voxW);
    if (uz >= VOXORDER) {
      uz = VOXORDER - 1;
    }

    AtomList* neighbourList = nullptr;
    for (int x = lx; x <= ux; x++) {
      for (int y = ly; y <= uy; y++) {
        for (int z = lz; z <= uz; z++) {
          for (AtomList* list = grid[x][y][z]; list; list = list->next) {
            for (unsigned int i = 0; i < list->count; i++) {
              const AtomRecord* temp = list->ptr[i];
              if (temp != &atom) {
                maxDist = range + temp->radius;
                if (within(atom, temp, maxDist)) {
                  insertAtomList(&neighbourList, temp);
                }
              }
            }
          }
        }
      }
    }

    return neighbourList;
  }

  bool testPoint(double* vect, double solvrad) {
    if (recordCache) {
      double dist = recordCache->radius + solvrad;
      double dx = recordCache->x - vect[0];
      double dy = recordCache->y - vect[1];
      double dz = recordCache->z - vect[2];
      if ((dx * dx + dy * dy + dz * dz) < (dist * dist)) {
        return false;
      }
      recordCache = nullptr;
    }

    for (const AtomList* list = neighbours; list; list = list->next) {
      for (unsigned int i = 0; i < list->count; i++) {
        const AtomRecord* ptr = list->ptr[i];
        double dist = ptr->radius + solvrad;
        double dx = ptr->x - vect[0];
        double dy = ptr->y - vect[1];
        double dz = ptr->z - vect[2];

        if ((dx * dx + dy * dy + dz * dz) < (dist * dist)) {
          recordCache = ptr;
          return false;
        }
      }
    }

    totalDots++;
    return true;
  }

  void determineCentreOfGravity(std::vector<AtomRecord>& memberAtoms) {
    double cx, cy, cz;

    cx = cy = cz = 0.0;
    for (const AtomRecord& atom : memberAtoms) {
      cx += atom.x;
      cy += atom.y;
      cz += atom.z;
    }

    unsigned int atomCount = memberAtoms.size();
    cenX = cx / atomCount;
    cenY = cy / atomCount;
    cenZ = cz / atomCount;
  }

  void generateStandardDots(int depth) {
    // Vertex co-ordinates of a unit icosahedron
    static const double Vertices[12][3] = {
        {0.00000000, -0.85065081, -0.52573111},
        {-0.52573111, 0.00000000, -0.85065081},
        {-0.85065081, -0.52573111, 0.00000000},
        {0.00000000, -0.85065081, 0.52573111},
        {0.52573111, 0.00000000, -0.85065081},
        {-0.85065081, 0.52573111, 0.00000000},
        {0.00000000, 0.85065081, -0.52573111},
        {-0.52573111, 0.00000000, 0.85065081},
        {0.85065081, -0.52573111, 0.00000000},
        {0.00000000, 0.85065081, 0.52573111},
        {0.52573111, 0.00000000, 0.85065081},
        {0.85065081, 0.52573111, 0.00000000}};

    // Face list of a unit icosahedron
    static const int Faces[20][4] = {
        {0, 1, 2},  {0, 1, 4},  {0, 2, 3},  {0, 3, 8},  {0, 4, 8},
        {1, 2, 5},  {1, 4, 6},  {1, 5, 6},  {2, 3, 7},  {2, 5, 7},
        {4, 8, 11}, {4, 6, 11}, {3, 8, 10}, {3, 7, 10}, {8, 10, 11},
        {5, 7, 9},  {5, 6, 9},  {6, 9, 11}, {7, 9, 10}, {9, 10, 11}};

    /* Count = 20*(4^Depth); */
    const unsigned long count = (20) * pow(4, depth);
    std::vector<DotStruct> dots(count);

    standardDots.radius = 1.0;
    standardDots.dots = dots;
    standardDots.count = 0;

    standardArea = 0.0;
    for (unsigned int i = 0; i < 20; i++)
      tesselate(&Vertices[Faces[i][0]][0], &Vertices[Faces[i][1]][0],
                &Vertices[Faces[i][2]][0], depth);
  }

  void generateSurfacePoints(int depth, bool typeFlag, double probeRadius,
                             int dotDensity) {
    if (!dotDensity) {
      generateStandardDots(depth);
    }

    if (typeFlag) {
      generateElemPoints(&elemA, 1.87, probeRadius, dotDensity);
      generateElemPoints(&elemC, 1.76, probeRadius, dotDensity);
      generateElemPoints(&elemN, 1.65, probeRadius, dotDensity);
      generateElemPoints(&elemO, 1.4, probeRadius, dotDensity);
      generateElemPoints(&elemX, 1.8, probeRadius, dotDensity);
    } else {
      generateElemPoints(&elemC, 1.87, probeRadius, dotDensity);
      generateElemPoints(&elemN, 1.5, probeRadius, dotDensity);
      generateElemPoints(&elemO, 1.4, probeRadius, dotDensity);
      generateElemPoints(&elemS, 1.84, probeRadius, dotDensity);
      generateElemPoints(&elemX, 1.44, probeRadius, dotDensity);
    }
  }

  void createVoxelGrid(int mask, std::vector<AtomRecord>& memberAtoms) {
    double minx, miny, minz;
    double maxx, maxy, maxz;

    minx = miny = minz = std::numeric_limits<double>::infinity();
    maxx = maxy = maxz = -std::numeric_limits<double>::infinity();

    freeList = nullptr;
    for (unsigned int x = 0; x < VOXORDER; x++) {
      for (unsigned int y = 0; y < VOXORDER; y++) {
        for (unsigned int z = 0; z < VOXORDER; z++) {
          grid[x][y][z] = nullptr;
        }
      }
    }
    // loop over atoms to find min + max x, y + z
    bool init = false;
    for (const AtomRecord& atom : memberAtoms) {
      if (!(atom.flag & mask)) {
        if (init) {
          if (atom.x > maxx) {
            maxx = atom.x;
          } else if (atom.x < minx) {
            minx = atom.x;
          }

          if (atom.y > maxy) {
            maxy = atom.y;
          } else if (atom.y < miny) {
            miny = atom.y;
          }

          if (atom.z > maxz) {
            maxz = atom.z;
          } else if (atom.z < minz) {
            minz = atom.z;
          }
        } else {
          maxx = minx = atom.x;
          maxy = miny = atom.y;
          maxz = minz = atom.z;
          init = true;
        }
      }
    }

    if (!init) return;

    voxX = VOXORDER / ((maxx - minx) + 0.1);
    voxU = minx;
    voxY = VOXORDER / ((maxy - miny) + 0.1);
    voxV = miny;
    voxZ = VOXORDER / ((maxz - minz) + 0.1);
    voxW = minz;

    for (AtomRecord& atom : memberAtoms) {
      if ((atom.flag & mask)) continue;

      // get grid positions and add to list
      unsigned int x = voxX * (atom.x - voxU);
      unsigned int y = voxY * (atom.y - voxV);
      unsigned int z = voxZ * (atom.z - voxW);
      insertAtomList(&grid[x][y][z], &atom);
    }
  }

  double surfaceArea(const AtomRecord& atom, const double solvrad,
                     const int dotDensity) {
    neighbours = findNeighbours(atom, atom.elem->radius + solvrad + solvrad);
    recordCache = nullptr;

    double surfacearea = 0.0;
    double vect[3];

    if (dotDensity) {
      for (const DotStruct& dot : atom.elem->dots) {
        vect[0] = atom.x + dot.v[0];
        vect[1] = atom.y + dot.v[1];
        vect[2] = atom.z + dot.v[2];
        if (testPoint(vect, solvrad)) {
          surfacearea += dot.area;
        }
      }
    } else {
      double factor = atom.elem->radius + solvrad;
      for (const DotStruct& dot : atom.elem->dots) {
        vect[0] = atom.x + factor * dot.v[0];
        vect[1] = atom.y + factor * dot.v[1];
        vect[2] = atom.z + factor * dot.v[2];
        if (testPoint(vect, solvrad)) {
          surfacearea += dot.area;
        }
      }
      surfacearea *= ((4.0 * M_PI) * atom.elem->radius2) / standardArea;
    }

    freeAtomList(neighbours);

    return surfacearea;
  }

  double calculateAccessibility(bool isProtein, double probeRadius,
                                int dotDensity,
                                std::vector<AtomRecord>& memberAtoms) {
    AtomRecord* ptr = nullptr;
    AtomRecord* prev = nullptr;
    unsigned int count = 0;
    double area = 0;

    totalDots = 0;
    double totalArea = 0.0;
    unsigned int atomCount = memberAtoms.size();

    bool init = false;

    if (isProtein) {
      for (unsigned int i = 0; i < atomCount; i++) {
        ptr = &memberAtoms[i];

        if (ptr->flag & (WaterFlag)) {
          ptr->bFactor = 0.0;
        } else if (!init) {
          if (!(ptr->flag & HetAtmFlag)) {
            init = true;
            prev = ptr;
            area = surfaceArea(*ptr, probeRadius, dotDensity);

            if (ptr->atmName != " OXT") {
              count = 1;
            } else {
              count = 0;
            }
            ptr->bFactor = (float)area;
          } else {
            ptr->bFactor = (float)0.0;
          }
        } else if (!sameResidue(ptr, prev)) {
          totalArea += area;
          checkResidue(prev, count);
          while (prev != ptr) {
            if (!(prev->flag & (HetAtmFlag))) {
              prev->bFactor = (float)area;
            }
            prev++;
          }

          if (!(ptr->flag & HetAtmFlag)) {
            area = surfaceArea(*ptr, probeRadius, dotDensity);
            if (ptr->atmName != " OXT") {
              count = 1;
            } else {
              count = 0;
            }
            ptr->bFactor = (float)area;
            prev = ptr;

          } else {
            ptr->bFactor = (float)0.0;
            init = false;
          }
        } else if (!(ptr->flag & HetAtmFlag)) {
          ptr->bFactor = (float)surfaceArea(*ptr, probeRadius, dotDensity);

          area += ptr->bFactor;
          if (ptr->atmName.compare(0, 4, " OXT") != 0) {
            count++;
          }
        }
      }

      if (init) {
        totalArea += area;
        checkResidue(prev, count);
        while (prev != ptr) {
          if (!(prev->flag & (HetAtmFlag))) {
            prev->bFactor = (float)area;
          }
          prev++;
        }
      }
    } else {
      for (const AtomRecord& atom : memberAtoms) {
        area = surfaceArea(atom, probeRadius, dotDensity);
        totalArea += area;
      }
    }
    return totalArea;
  };

  double partialVolume(const AtomRecord& atom, const double solvrad) {
    double rad = atom.elem->radius + solvrad;
    neighbours = findNeighbours(atom, rad + solvrad);
    recordCache = nullptr;

    unsigned int count = 0;

    double px, py, pz;
    px = py = pz = 0.0;

    double vect[3];
    for (const DotStruct& dot : atom.elem->dots) {
      vect[0] = atom.x + rad * dot.v[0];
      vect[1] = atom.y + rad * dot.v[1];
      vect[2] = atom.z + rad * dot.v[2];
      if (testPoint(vect, solvrad)) {
        px += dot.v[0];
        py += dot.v[1];
        pz += dot.v[2];
        count++;
      }
    }
    freeAtomList(neighbours);

    /* Calculate Volume with Gauss-Ostrogradskii Theorem */
    double partialvol =
        (atom.x - cenX) * px + (atom.y - cenY) * py + (atom.z - cenZ) * pz;
    partialvol = rad * rad * (partialvol + rad * count);

    return partialvol;
  }

  double calculateVolume(const double probeRadius,
                         const std::vector<AtomRecord>& memberAtoms) {
    double totalvol = 0.0;
    for (const AtomRecord& atom : memberAtoms) {
      totalvol += partialVolume(atom, probeRadius);
    }
    totalvol *= ((4.0 / 3) * M_PI) / standardDots.count;

    return totalvol;
  };
};

// constructor definition
DoubleCubicLatticeVolume::DoubleCubicLatticeVolume(const ROMol* mol,
                                                   bool isProtein,
                                                   bool includeLigand,
                                                   double probeRadius,
                                                   int depth, int dotDensity) {
  //! Class for calculation of the Shrake and Rupley surface area and volume
  //! using the Double Cubic Lattice Method.
  //!
  //! Frank Eisenhaber, Philip Lijnzaad, Patrick Argos, Chris Sander and
  //! Michael Scharf, "The Double Cubic Lattice Method: Efficient Approaches
  //! to Numerical Integration of Surface Area and Volume and to Dot Surface
  //! Contouring of Molecular Assemblies", Journal of Computational Chemistry,
  //! Vol. 16, No. 3, pp. 273-284, 1995.

  /*!

    \param mol: input molecule or protein
    \param includeLigand: flag to trigger inclusion of bound
    ligand in surface area and volume calculations where molecule is a protein
    [default false]
    \param probeRadius: radius of the sphere representing the probe solvent
    atom
    \param depth: controls the number of dots per atom
    \param dotDensity: controls density of dots per atom
    \return class
    object
  */

  if (depth > maxDepth) {
    throw std::range_error("Error: supplied depth exceeds maximum");
  }

  if (dotDensity > maxDensity) {
    throw std::range_error("Error: supplied density exceeds maximum");
  }

  if (depth < 0) {
    depth = 0;
  }

  if (probeRadius < 0.0) {
    probeRadius = (isProtein) ? 1.4 : 1.2;
  }

  // default 320 faces, protein
  if (!isProtein && probeRadius == 1.4 && depth == 2) {
    probeRadius = 1.2;
    depth = 4;  // 5120 faces, ligand
  }

  mol = MolOps::removeAllHs(*mol, false);

  State s;
  s.generateSurfacePoints(depth, isProtein, probeRadius, dotDensity);

  const Conformer& conf = mol->getConformer();

  std::vector<AtomRecord> memberAtoms;

  for (const Atom* atom : mol->atoms()) {
    AtomRecord curr_atom(atom, conf);
    curr_atom.initFlag();
    if (!curr_atom.solventFlag) {
      curr_atom.elem = s.getAtomElem(curr_atom.atmName, isProtein);
      curr_atom.radius = curr_atom.elem->radius;
      memberAtoms.push_back(curr_atom);
    }
  }

  int mask = includeLigand ? WaterFlag : HetAtmFlag;
  s.determineCentreOfGravity(memberAtoms);
  s.createVoxelGrid(mask, memberAtoms);
  surfaceArea =
      s.calculateAccessibility(isProtein, probeRadius, dotDensity, memberAtoms);
  totalVolume = s.calculateVolume(probeRadius, memberAtoms);
  vdwVolume = s.calculateVolume(0.0, memberAtoms);
  compactness = calculateCompactness(surfaceArea, totalVolume);
  packingDensity = (vdwVolume / totalVolume);
}

}  // namespace Descriptors
}  // namespace RDKit