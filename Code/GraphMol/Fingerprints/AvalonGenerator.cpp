//
//  Copyright (C) 2025 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  This is a pure RDKit port of the fingerprinting algorithm from the Avalon
//  toolkit (SetFingerprintBits() / SetFingerprintCountsWithFocus() in
//  ssmatch.c). The feature enumeration and hashing follow the original code
//  as closely as possible, but molecule perception (aromaticity, rings,
//  hydrogen counts) is done with RDKit's own machinery.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Fingerprints/AvalonGenerator.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <RDGeneral/Invariant.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace RDKit {
namespace AvalonFP {
namespace {

// ---- minimal stand-ins for the data structures used by the Avalon code ----
#define TRUE 1
#define FALSE 0
#define NONE 0
#define ZERO_COUNT 1
#define SINGLE 1
#define DOUBLE 2
#define TRIPLE 3
#define AROMATIC 4
#define ANY_BOND 8
#define SUB_ONE 1
#define SUB_MORE 6
#define SUB_AS_IS -2

#define USE_RING_PATTERN 0x000001
#define USE_RING_PATH 0x000002
#define USE_ATOM_SYMBOL_PATH 0x000004
#define USE_ATOM_CLASS_PATH 0x000008
#define USE_ATOM_COUNT 0x000010
#define USE_AUGMENTED_ATOM 0x000020
#define USE_HCOUNT_PATH 0x000040
#define USE_HCOUNT_CLASS_PATH 0x000080
#define USE_HCOUNT_PAIR 0x000100
#define USE_BOND_PATH 0x000200
#define USE_AUGMENTED_BOND 0x000400
#define USE_RING_SIZE_COUNTS 0x000800
#define USE_DEGREE_PATH 0x001000
#define USE_CLASS_SPIDERS 0x002000
#define USE_FEATURE_PAIRS 0x004000
#define USE_SCAFFOLD_IDS 0x100000
#define USE_SCAFFOLD_COLORS 0x200000
#define USE_SCAFFOLD_LINKS 0x400000
#define USE_NON_SSS_BITS 0xF00000

struct reaccs_atom_t {
  std::string atom_symbol;
  int color = 0;
  int rsize_flags = 0;
  int sub_desc = NONE;
  int charge = 0;
};

struct reaccs_bond_t {
  std::array<int, 2> atoms{};  // 1-based atom numbers
  int bond_type = 0;
  int color = 0;
  int rsize_flags = 0;
};

struct reaccs_molecule_t {
  int n_atoms = 0;
  int n_bonds = 0;
  reaccs_atom_t *atom_array = nullptr;
  reaccs_bond_t *bond_array = nullptr;
};

struct neighbourhood_t {
  int n_ligands = 0;
  std::vector<int> atoms;  // 0-based atom indices
  std::vector<int> bonds;  // 0-based bond indices
};

// ---- hashing ----
uint64_t next_hash(uint64_t hash, uint64_t data) {
  hash += data;
  hash += (uint64_t)(hash << (uint64_t)10);
  hash ^= (uint64_t)(hash >> (uint64_t)6);
  return hash;
}

uint64_t hash_position(uint64_t hash, int nslots) {
  hash += (uint64_t)(hash << (uint64_t)3);
  hash ^= (uint64_t)(hash >> (uint64_t)11);
  hash += (uint64_t)(hash << (uint64_t)15);
  return (hash % (uint64_t)nslots);
}

int AtomicNumberFromSymbol(const std::string &symbol) {
  if (symbol == "*") {
    return 0;
  }
  return PeriodicTable::getTable()->getAtomicNumber(symbol.c_str());
}

// true if symbol is one of the comma separated tokens in list
int AtomSymbolMatch(std::string_view symbol, std::string_view list) {
  while (!list.empty()) {
    const auto comma = list.find(',');
    if (list.substr(0, comma) == symbol) {
      return TRUE;
    }
    if (comma == std::string_view::npos) {
      break;
    }
    list.remove_prefix(comma + 1);
  }
  return FALSE;
}

// ---- the Avalon algorithm ----
#define RING_PATTERN_SEED 11
#define RING_PATH_SEED 13
#define ATOM_SYMBOL_PATH_SEED 17
#define ATOM_CLASS_PATH_SEED 23
#define ATOM_COUNT_SEED 31
#define AUGMENTED_ATOM_SEED 37
#define HCOUNT_PATH_SEED 41
#define HCOUNT_CLASS_PATH_SEED 43
#define HCOUNT_PAIR_SEED 47
#define BOND_PATH_SEED 53
#define AUGMENTED_BOND_SEED 61
#define RING_SIZE_SEED 67
#define DEGREE_PATH_SEED 71
#define CLASS_SPIDER_SEED 79
#define RING_CLOSURE_SEED 101
#define NON_SSS_SEED 179
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* new macro to convert the current seed value into the 'incremented' one */
#define NEXT_SEED(seed, increment) next_hash(seed, increment)
#define ADD_BIT(counts, ncounts, seed) (counts[hash_position(seed, ncounts)]++)
#define ADD_BIT_COUNT(counts, ncounts, seed, count) \
  (counts[hash_position(seed, ncounts)] += count)
#define SET_BIT(bytes, nbytes, seed) \
  (bytes[((seed) / 8) % nbytes] |= 0xFF & (1 << ((seed) % 8)))

/* Flags to be used to control recursive processing */
#define PROCESS_RING_CLOSURES 0x0001
#define PROCESS_CHAINS 0x0002
#define FORCED_HETERO_END 0x0004
#define IGNORE_PATH_SYMBOL 0x0008
#define IGNORE_TERM_SYMBOL 0x0010
#define FORCED_RING_PATH 0x0020
#define STOP_AT_HEAVY_ATOM 0x0040
#define DEBUG_PATH 0x0100

#define ANY_COLOR 113

#define CSP3 19
#define HETERO 23
#define GENERIC (-1)

#define SPECIAL_RING (0xFC & ~(1 << 6))

static void SetPathLengthFlags(struct reaccs_molecule_t *mp,
                               std::vector<int> &touched_indices, int start_index,
                               int path_length, int current_index, int max_size,
                               std::vector<std::vector<int>> &length_matrix,
                               const std::vector<neighbourhood_t> &nbp,
                               int exclude_atom)
/*
 * Recursively traces the neighbouring of an atom (start_index+1)
 * collecting the path_lengths as bit flags in length_matrix[][].
 */
{
  for (int i = 0; i < nbp[current_index].n_ligands; i++) {
    if (path_length + 1 > max_size) { /* don't go too far */
      continue;
    }
    const auto ai = nbp[current_index].atoms[i];
    if (ai + 1 == exclude_atom) {
      continue;
    }
    if (touched_indices[ai]) {
      continue; /* don't walk backwards */
    }
    if (mp->atom_array[ai].color == 0) {
      continue;
    }
    touched_indices[ai] = 1; /* updating */
    length_matrix[start_index][ai] |= 1 << (path_length + 1);
    SetPathLengthFlags(mp, touched_indices, start_index, path_length + 1, ai,
                       max_size, length_matrix, nbp, exclude_atom);
    touched_indices[ai] = 0; /* down-dating */
  }
}

static void SpecialNeighboursRec(
    struct reaccs_molecule_t *mp, std::vector<int> &touched_indices,
    int path_length,
    int current_index, int max_size,
    /* count of sp3 carbons with >= 3 C neighbours */
    int csp3[],
    /* count of hetero atoms */
    int hetero[], const std::vector<neighbourhood_t> &nbp, int exclude_atom)
/*
 * Recursively traces the neighbouring of an atom collecting
 * counts of special atoms at certain graph distances.
 *
 * exclude_atom terminates neighbourhood search paths.
 */
{
  for (int i = 0; i < nbp[current_index].n_ligands; i++) {
    if (path_length + 1 > max_size) { /* don't go too far */
      continue;
    }
    const auto ai = nbp[current_index].atoms[i];
    if (ai + 1 == exclude_atom) {
      continue;
    }
    if (touched_indices[ai]) {
      continue; /* don't walk backwards */
    }
    if (mp->atom_array[ai].color <= 0) {
      continue;
    }
    if (mp->atom_array[ai].color == CSP3) {
      csp3[path_length]++;
    } else if (mp->atom_array[ai].color == HETERO) {
      hetero[path_length]++;
    }
    /* only carbon-connected SPIDERS (beta atoms are always included) */
    if (mp->atom_array[ai].color == HETERO && path_length > 1) {
      continue;
    }
    touched_indices[ai] = 1; /* updating */
    SpecialNeighboursRec(mp, touched_indices, path_length + 1, ai, max_size,
                         csp3, hetero, nbp, exclude_atom);
    touched_indices[ai] = 0; /* down-dating */
  }
}
int SetPathBitsRec(struct reaccs_molecule_t *mp,
                   const std::vector<neighbourhood_t> &nbp,
                   int *fp_counts, int ncounts, uint64_t seed,
                   std::vector<int> &touched_indices, int nbonds,
                   int minbonds, int maxbonds, int sprout_index,
                   int first_index, int last_index, int flags, int exclude_atom)
/*
 * Recursively enumerates the paths through *mp. The next sprouting
 * step is done on the atom (sprout_index+1). seed represents the
 * hash_code processing value up to and including this atom. touched_indices[i]
 * is > 0 if atom (i+1) is already included in the parent path.
 *
 * Paths through exclude_atom are terminated.
 */
{
  int result;
  int ai, bi, acolor, bcolor;
  uint64_t old_seed;
  struct reaccs_atom_t *ap;

  result = 0;
  old_seed = seed;
  if (nbonds > maxbonds) {
    return (result);
  }
  for (int i = 0; i < nbp[sprout_index].n_ligands; i++) {
    ai = nbp[sprout_index].atoms[i];
    if (ai == last_index) {
      continue;
    }
    if (ai + 1 == exclude_atom) {
      continue;
    }
    bi = nbp[sprout_index].bonds[i];
    ap = mp->atom_array + ai;
    if ((flags & FORCED_RING_PATH) && mp->bond_array[bi].rsize_flags == 0) {
      continue;
    }
    if (ap->color >= 18 && ap->color != ANY_COLOR &&
        0 != (flags & STOP_AT_HEAVY_ATOM)) {
      continue;
    }
    if (touched_indices[ai] > 0) /* ring closure */
    {
      if (0 == (flags & PROCESS_RING_CLOSURES)) {
        continue;
      }
      if (touched_indices[ai] > 1) {
        continue;  // not just a plain ring
      }
      bcolor = mp->bond_array[bi].color;
      if (bcolor == 0) {
        continue;
      }
      acolor = mp->atom_array[ai].color;
      if (acolor == 0 && 0 == (flags & IGNORE_PATH_SYMBOL)) {
        continue;
      }
      old_seed = seed;
      seed = NEXT_SEED(seed, bcolor * 16);
      if (0 == (flags & IGNORE_PATH_SYMBOL)) {
        seed = NEXT_SEED(seed, acolor);
      }
      seed = NEXT_SEED(
          seed, touched_indices[ai]);  // codes how far this closure points back
      ADD_BIT(fp_counts, ncounts, seed);
      seed = NEXT_SEED(seed,
                       RING_CLOSURE_SEED *
                           (nbonds - touched_indices[ai]));  // codes ring size
      ADD_BIT(fp_counts, ncounts, seed);
      result++;
      seed = old_seed;
    } else /* normal path */
    {
      bcolor = mp->bond_array[bi].color;
      if (bcolor == 0) {
        continue;
      }
      acolor = mp->atom_array[ai].color;
      if (acolor == 0 && 0 == (flags & IGNORE_PATH_SYMBOL)) {
        continue;
      }
      /* saving and updating */
      touched_indices[ai] = nbonds + 1;
      old_seed = seed;
      seed = NEXT_SEED(seed, bcolor * 16);
      if (0 == (flags & IGNORE_PATH_SYMBOL)) {
        seed = NEXT_SEED(seed, acolor);
      }
      if (nbonds >= minbonds &&
          (acolor > 0 || 0 != (flags & IGNORE_TERM_SYMBOL)) &&
          0 != (flags & PROCESS_CHAINS)) {
        if (0 == (flags & FORCED_HETERO_END) ||
            0 != (flags & IGNORE_TERM_SYMBOL) || acolor != 6) {
          if (acolor > 1) /* don't hit hydrogens! */
          {
            if (0 != (flags & IGNORE_TERM_SYMBOL)) {
              ADD_BIT(fp_counts, ncounts, seed);
            } else {
              ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 17 * acolor));
            }
            result++;
          }
        }
      }
      if (nbonds + 1 <= maxbonds) /* continue recursion */
      {
        result +=
            SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices,
                           nbonds + 1, minbonds, maxbonds, ai, first_index,
                           sprout_index, flags, exclude_atom);
      }

      /* restoring and down dating */
      touched_indices[ai] = 0;
      seed = old_seed;
    }
  }

  return result;
}

#define HETERO_FLAG 0x0100
#define RING_SUBST_FLAG 0x0200
#define QUART_FLAG 0x0400
#define CSP3_FLAG 0x0800
#define RS_SPECIAL_FLAG 0x1000
#define TYPE_MASK 0x00FF
#define C_FLAG 0x0001
#define O_FLAG 0x0002
#define N_FLAG 0x0003
#define S_FLAG 0x0004
#define P_FLAG 0x0005
#define X_FLAG 0x0006

int SetFeatureBits(struct reaccs_molecule_t *mp, int *fp_counts, int ncounts,
                   int start_flags, int end_flags, int path_min, int path_max,
                   int use_counts, int use_atom_types,
                   const std::vector<std::vector<int>> &length_matrix,
                   uint64_t start_seed, int exclude_atom) {
  int result = 0;
  int coli, colj;
  uint64_t seed_i, seed;
  std::vector<int> counts(ncounts * 4, 0);
  for (int i = 0; i < mp->n_atoms; i++) {
    if (i + 1 == exclude_atom) {
      continue;
    }
    coli = mp->atom_array[i].color;
    if (0 == (coli & start_flags)) {
      continue;
    }
    if (use_atom_types) {
      if (0 == (coli & TYPE_MASK)) {
        continue;  // ignore generic atoms
      }
      seed_i = NEXT_SEED(start_seed, coli & TYPE_MASK);
    } else {
      seed_i = start_seed;
    }
    for (int j = 0; j < mp->n_atoms; j++) {
      if (j + 1 == exclude_atom) {
        continue;
      }
      colj = mp->atom_array[j].color;
      if (0 == (colj & end_flags)) {
        continue;
      }
      if (use_atom_types) {
        if (0 == (colj & TYPE_MASK)) {
          continue;  // ignore generic atoms
        }
        seed = NEXT_SEED(seed_i, colj & TYPE_MASK);
      } else {
        seed = seed_i;
      }
      for (int k = path_min; k <= path_max; k++) {
        if ((1 << k) & length_matrix[i][j]) {
          /* count the features */
          counts[(k * 19 + seed) % (ncounts * 4)]++;
          // ADD_BIT(fp_counts, ncounts, k*19+seed);
          result++;
        }
      }
    }
  }
  /* Set the bits */
  for (int i = 0; i < ncounts * 4; i++) {
    if (counts[i] > 0) {
      ADD_BIT(fp_counts, ncounts, i);
      result++;
    }
    if (use_counts && counts[i] > 1) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i, 19), counts[i]);
      result++;
    }
    if (use_counts && counts[i] > 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i, 29), counts[i]);
      result++;
    }
    if (use_counts && counts[i] > 4) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i, 59), counts[i]);
      result++;
    }
    if (use_counts && counts[i] > 8) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i, 79), counts[i]);
      result++;
    }
  }
  return (result);
}
int CountFingerprintPatterns(
    reaccs_molecule_t *mp, const std::vector<neighbourhood_t> &nbp,
    int *H_count, int *atom_status, int *bond_status, int *fp_counts,
    int ncounts, int which_bits, int as_query, int exclude_atom) {
  uint64_t seed, old_seed;
  std::vector<int> touched_indices(mp->n_atoms, 0);
  std::vector<int> degree(mp->n_atoms, 0);
  std::vector<int> cdegree(mp->n_atoms, 0);
  std::vector<int> unsaturated(mp->n_atoms, 0);
  std::vector<int> nspecial(mp->n_atoms, 0);
  std::vector<int> extcon;
  std::vector<int> extcon2;
  std::vector<std::vector<int>> length_matrix;
  reaccs_atom_t *ap, *ap1, *ap2, *ap3;
  reaccs_bond_t *bp;
  int tmp;

  uint64_t prod, sum, sumi, prodi, sumj, prodj;
  int ai, ai1, ai2;
  int qq_count;
  int qc_count;
  int nrbonds;
  int nrare_atoms;
  int is_rare;
  int nqtmp, nmulti;
  int hash;
  int nbits;
  int result; /* the number of paths enumerated */
  constexpr int NCOUNT_HASH = 128;
  constexpr int NCOUNT_SEED_HASH = 128 * 128;
  std::array<int, NCOUNT_HASH> atom_type_count_hash{};
  std::array<int, NCOUNT_SEED_HASH> atom_type_count_seed_hash{};
  constexpr int MAX_SPIDER = 7;
  std::array<int, MAX_SPIDER + 1> csp3{};
  std::array<int, MAX_SPIDER + 1> hetero{};
  int tmp1, tmp2;
  int ndouble, naromatic;
  std::array<std::array<int, 15>, 15> rscounts{};
  int nringch2, nfusionch, nspiro, nfusionb;
  int flags;
  int changed;

  result = 0;
  nrare_atoms = 0;
  /* Set the color property to represent all different atom types */
  ap = mp->atom_array;
  for (int i = 0; i < mp->n_atoms; i++, ap++) {
    unsaturated[i] = FALSE;
    ap->color = AtomicNumberFromSymbol(ap->atom_symbol);
    if (ap->color <= 1) {
      ap->color = 0; /* ignore hydrogens */
    }
    /* mark special atom types */
    if (ap->color > 115) {
      ap->color = -1;
    }
    if (ap->atom_symbol == "A") {
      ap->color = -1;
    }
    is_rare =
        ap->color > 0 && !AtomSymbolMatch(ap->atom_symbol, "C,H,O,N,S,P,Cl,F");
    if (is_rare) {
      if (exclude_atom != i + 1 || exclude_atom <= 0) {
        nrare_atoms++;
      }
    }
  }

  ndouble = 0;
  naromatic = 0;
  nfusionb = 0;
  /* Set the color property to represent the different bond type classes */
  bp = mp->bond_array;
  for (int i = 0; i < mp->n_bonds; i++, bp++) {
    if (bp->bond_type == SINGLE) {
      bp->color = 1;
    } else if (bp->bond_type == DOUBLE) {
      bp->color = 2;
      if (bp->atoms[0] != exclude_atom && bp->atoms[1] != exclude_atom) {
        ndouble++;
      }
    } else if (bp->bond_type == TRIPLE) {
      bp->color = 3;
    } else if (bp->bond_type == AROMATIC) {
      bp->color = 4;
      if (bp->atoms[0] != exclude_atom && bp->atoms[1] != exclude_atom) {
        naromatic++;
      }
    } else {
      bp->color = 0;
    }
    if (bp->color > 1) {
      unsaturated[bp->atoms[0] - 1] = TRUE;
      unsaturated[bp->atoms[1] - 1] = TRUE;
    }

    /* Count non-hydrogen degree */
    if (mp->atom_array[bp->atoms[0] - 1].color != 0 &&
        mp->atom_array[bp->atoms[1] - 1].color != 0) {
      degree[bp->atoms[0] - 1]++;
      degree[bp->atoms[1] - 1]++;
    }
    /* Count carbon degree */
    if (bp->bond_type == DOUBLE) {
      nspecial[bp->atoms[0] - 1]++;
      nspecial[bp->atoms[1] - 1]++;
    } else if (bp->bond_type == TRIPLE) {
      nspecial[bp->atoms[0] - 1] += 2;
      nspecial[bp->atoms[1] - 1] += 2;
    }
    if (bp->bond_type != SINGLE) {
      continue;
    }
    if (mp->atom_array[bp->atoms[0] - 1].color != 0 &&
        mp->atom_array[bp->atoms[1] - 1].color == 6) {
      cdegree[bp->atoms[0] - 1]++;
    }
    if (mp->atom_array[bp->atoms[1] - 1].color != 0 &&
        mp->atom_array[bp->atoms[0] - 1].color == 6) {
      cdegree[bp->atoms[1] - 1]++;
    }

    if (exclude_atom <= 0 ||
        (bp->atoms[0] != exclude_atom && bp->atoms[1] != exclude_atom)) {
      if (atom_status[bp->atoms[0] - 1] > 2 &&
          atom_status[bp->atoms[1] - 1] > 2) {  // ring fusion
        nfusionb++;
      }
    }
  }
  /* ignore special atom types for further processing */
  ap = mp->atom_array;
  for (int i = 0; i < mp->n_atoms; i++, ap++) {
    if (ap->color < 0) {
      ap->color = 0;
    }
  }

  if (which_bits & USE_ATOM_COUNT) {
    /* Collect hashed counts of atom types with hydrogen counts */
    std::fill(atom_type_count_hash.begin(), atom_type_count_hash.end(), 0);
    std::fill(atom_type_count_seed_hash.begin(),
              atom_type_count_seed_hash.end(), 0);
    nringch2 = 0;
    nfusionch = 0;
    nspiro = 0;
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (ap->color == 0) {
        continue;
      }
      if (ap->color == 6) {
        /* count CH2 in ring */
        if (atom_status[i] > 0 && H_count[i + 1] >= 2) {
          nringch2++;
        }
        /* count ring fusion CH atoms */
        if (atom_status[i] > 2 && H_count[i + 1] >= 1) {
          nfusionch++;
        }
        if (atom_status[i] > 3) {
          nspiro++;
        }
        continue;  // carbon has a special count processing
      }
      hash = 0;
      seed = NEXT_SEED(317 * ATOM_COUNT_SEED, 507);
      seed = NEXT_SEED(seed, ap->color + 17);
      atom_type_count_seed_hash[seed % NCOUNT_SEED_HASH]++;
      /* normal hetero only with hydrogen */
      if ((ap->color == 7 || ap->color == 8) && H_count[i + 1] <= 0) {
        continue;
      }
      hash = hash * 7 + ap->color + 13;
      atom_type_count_hash[hash % NCOUNT_HASH]++;
      seed = NEXT_SEED(seed, ap->color + 13);
      atom_type_count_seed_hash[seed % NCOUNT_SEED_HASH]++;
      /* one more bit for rare types */
      if (ap->color != 7 && ap->color != 8) {
        hash = hash * 7 + ap->color + 2 * 13;
        atom_type_count_hash[hash % NCOUNT_HASH]++;
        seed = NEXT_SEED(seed, ap->color + 2 * 13);
        atom_type_count_seed_hash[seed % NCOUNT_SEED_HASH]++;
      }
    }
    nbits = 0;
    /* Now, we set the corresponding bits */

    if (ncounts <= 2048) {  // old atom count fingerprints
      for (int i = 0; i < NCOUNT_HASH; i++) {
        if (atom_type_count_hash[i] > 0) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 19, 3),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 1) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 23, 5),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 2) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 29, 7),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 4) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 31, 11),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 8) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 37, 13),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 16) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 41, 17),
                        atom_type_count_hash[i]);
          nbits++;
        }
      }
    } else {  // new atom count fingerprints
      for (int i = 0; i < NCOUNT_SEED_HASH; i++) {
        if (atom_type_count_seed_hash[i] > 0) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 19, 3),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 1) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 23, 5),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 2) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 29, 7),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 3) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 29, 71),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 4) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 31, 11),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 6) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 31, 113),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 8) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 37, 13),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 12) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 37, 133),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 16) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 41, 17),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
      }
    }

    if (naromatic >= 10) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 341),
                    naromatic);
      nbits++;
    }
    if (naromatic >= 14) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 441),
                    naromatic);
      nbits++;
    }
    if (naromatic >= 18) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 541),
                    naromatic);
      nbits++;
    }
    if (naromatic >= 22) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 641),
                    naromatic);
      nbits++;
    }
    if (naromatic >= 26) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 741),
                    naromatic);
      nbits++;
    }

    if (ndouble >= 3) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 371),
                    naromatic);
      nbits++;
    }
    if (ndouble >= 5) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 471),
                    naromatic);
      nbits++;
    }
    if (ndouble >= 8) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 571),
                    naromatic);
      nbits++;
    }

    if (nringch2 >= 6) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 411),
                    nringch2);
      nbits++;
    }
    if (nringch2 >= 12) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 511),
                    nringch2);
      nbits++;
    }
    if (nringch2 >= 22) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 611),
                    nringch2);
      nbits++;
    }

    if (nfusionch >= 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 421),
                    nfusionch);
      nbits++;
    }
    if (nfusionch >= 4) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 521),
                    nfusionch);
      nbits++;
    }

    // Note spiro atoms
    if (nspiro >= 1) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 322),
                    nspiro);
      nbits++;
    }
    if (nspiro >= 1) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 422),
                    nspiro);
      nbits++;
    }
    if (nspiro >= 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 522),
                    nspiro);
      nbits++;
    }
    if (nspiro >= 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 622),
                    nspiro);
      nbits++;
    }

    if (nfusionb >= 1) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 425),
                    nfusionb);
      nbits++;
    }
    if (nfusionb >= 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 525),
                    nfusionb);
      nbits++;
    }
    if (nfusionb >= 3) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 625),
                    nfusionb);
      nbits++;
    }
    if (nfusionb >= 5) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 825),
                    nfusionb);
      nbits++;
    }
    if (nrare_atoms > 0) {
      seed = NEXT_SEED(317 * ATOM_COUNT_SEED, 5);
      ADD_BIT_COUNT(fp_counts, ncounts, seed, nrare_atoms);
      nbits++;
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed, 101), nrare_atoms);
      nbits++;
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed, 103), nrare_atoms);
      nbits++;
      if (nrare_atoms > 1) {
        ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed, 203), nrare_atoms);
        nbits++;
      }
      if (nrare_atoms > 3) {
        ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed, 211), nrare_atoms);
        nbits++;
      }
    }
    result += nbits;
  }

  if (which_bits & USE_ATOM_SYMBOL_PATH) {
    seed = ATOM_SYMBOL_PATH_SEED;
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (ap->color <= 0) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      /* Ignore common atom types */
      if (ap->color != 6 && ap->color != 7 && ap->color != 8) {
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
      }
      /* fingerprint two atom pairs if substitution count is defined */
      if (1) {  // [TODO] 1
        if (ap->color != 6 &&
            (!as_query || ap->sub_desc == SUB_AS_IS ||
             (ap->sub_desc != NONE && ap->sub_desc != SUB_MORE &&
              ap->sub_desc == degree[i] + SUB_ONE - 1))) {
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed + 12347, touched_indices, 1, 1,
              1, /* path length 1 to 1 */
              i, 0, -1, FORCED_HETERO_END | PROCESS_CHAINS, exclude_atom);
        }
      }
      /* 2 bond paths for not very common starts */
      if (1) {  // [TODO] 2
        if (ap->color >= 10 || (ap->color == 7 && degree[i] > 0) ||
            (ap->color == 8 && degree[i] > 1) ||
            (ap->color == 6 && degree[i] > 2 && atom_status[i] > 0)) {
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 1,
              2, /* path length 1 to 2 */
              i, 0, -1,
              STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  FORCED_HETERO_END | PROCESS_CHAINS,
              exclude_atom);
        }
      }
      /* Add more paths starting at special atoms */
      if (1) {  // [TODO] 3
        if (ap->color > 6) {
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, 217 * ap->color + seed,
              touched_indices, 1, 3, 4, /* path length 3 to 4 */
              i, 0, -1,
              IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES |
                  STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  PROCESS_CHAINS,
              exclude_atom);

          if (ap->color > 10 &&
              ap->color <= 18) {  // only third row of periodic table
            result += SetPathBitsRec(
                mp, nbp, fp_counts, ncounts, 17 + seed, touched_indices, 1, 5,
                7, i, 0, -1,
                FORCED_HETERO_END | IGNORE_PATH_SYMBOL |
                    STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                    // PROCESS_RING_CLOSURES |      // rings containing metal
                    // cause too many bits
                    PROCESS_CHAINS,
                exclude_atom);
          }
        }
      }
      if (1) {  // [TODO] 4
        if ((ap->color == 7 ||
             ap->color == 8) &&  // only do this for common hetero elements
            cdegree[i] > 2) {
          seed = NEXT_SEED(seed, ap->color * 23);
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 2,
              6, /* path length 1 to 4 */
              i, 0, -1,
              STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  IGNORE_TERM_SYMBOL |
                  // IGNORE_PATH_SYMBOL |
                  PROCESS_CHAINS,
              exclude_atom);
        }
      }

      seed = old_seed;
      /* Add another set of bits for paths starting at more crowded atoms */
      if (1) {                                   // [TODO] 5
        if (degree[i] >= 4 && ap->color != 5 &&  // don't do it for boron
            ap->color < 18)  // don't do it for transition metals
        {
          seed = NEXT_SEED(seed, ap->color);
          if (1) {
            result += SetPathBitsRec(
                mp, nbp, fp_counts, ncounts, NEXT_SEED(seed, 3 * 107),
                touched_indices, 1, 2, 4, /* path length 2 to 4 */
                i, 0, -1,
                STOP_AT_HEAVY_ATOM |     // don't cross very heavy atoms
                    IGNORE_TERM_SYMBOL | /* gedeck.mol */
                    PROCESS_CHAINS,
                exclude_atom);
          }
        }
      }

      seed = old_seed;
      /* Add another set of bits for paths starting at spiro atoms */
      if (1) {                                        // [TODO] 6
        if (atom_status[i] >= 4 && ap->color != 5 &&  // don't do it for boron
            ap->color < 18)  // don't do it for transition metals
        {
          seed = NEXT_SEED(seed, ap->color + 55);
          // if (1)
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, NEXT_SEED(seed, 3 * 109),
              touched_indices, 1, 2, 4, /* path length 1 to 4 */
              i, 0, -1,
              STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  IGNORE_TERM_SYMBOL | IGNORE_PATH_SYMBOL |
                  PROCESS_RING_CLOSURES | PROCESS_CHAINS,
              exclude_atom);
        }
      }

      seed = old_seed;
      /* Add bits for paths starting with rare bond orders */
      if (1) {
        for (int j = 0; j < nbp[i].n_ligands; j++) {
          bp = &mp->bond_array[nbp[i].bonds[j]];
          ai = nbp[i].atoms[j];
          if (ai + 1 == exclude_atom) {
            continue;
          }
          if (bp->color == 0) {
            continue;
          }
          if (bp->bond_type != DOUBLE && bp->bond_type != TRIPLE) {
            continue;
          }
          seed = old_seed;
          seed = NEXT_SEED(seed, ap->color);
          seed = NEXT_SEED(seed, bp->color * 613);
          touched_indices[ai] = 1; /* updating */
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 2, 5, 5, ai,
              0, i,
              STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  IGNORE_PATH_SYMBOL |
                  // IGNORE_TERM_SYMBOL |
                  PROCESS_RING_CLOSURES | (0 * PROCESS_CHAINS),
              exclude_atom);
          touched_indices[ai] = 0; /* down-dating */
        }
      }

      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
  }

  if (which_bits & USE_AUGMENTED_ATOM) {
    /* Set bits for all triples of atoms connected to a common atom */
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }

      /* Set bit for atoms with more than one double or a triple bond */
      if (nspecial[i] >= 2) {
        ADD_BIT(fp_counts, ncounts,
                NEXT_SEED(ap->color * AUGMENTED_ATOM_SEED, 101));
        // ADD_BIT(fp_counts, ncounts,
        // NEXT_SEED(ap->color*AUGMENTED_ATOM_SEED,301));
        result += 1;
      }
      old_seed = seed;
      // Add some bits for hydrogen counted or hetero central atoms with
      // hetero neighbours
      if ((H_count[i + 1] > 0 || ap->color != 6) && degree[i] >= 2) {
        for (int i1 = 0; i1 < nbp[i].n_ligands; i1++) {
          if (mp->atom_array[nbp[i].atoms[i1]].color == 0) {
            continue;
          }
          if (nbp[i].atoms[i1] + 1 == exclude_atom) {
            continue;
          }
          if (mp->bond_array[nbp[i].bonds[i1]].color == 0) {
            continue;
          }
          for (int i2 = i1 + 1; i2 < nbp[i].n_ligands; i2++) {
            seed = NEXT_SEED(AUGMENTED_ATOM_SEED, 97);
            seed = NEXT_SEED(seed, ap->color);
            if (mp->atom_array[nbp[i].atoms[i2]].color == 0) {
              continue;
            }
            if (nbp[i].atoms[i2] + 1 == exclude_atom) {
              continue;
            }
            if (mp->bond_array[nbp[i].bonds[i2]].color == 0) {
              continue;
            }
            if (mp->atom_array[nbp[i].atoms[i1]].color == 6 &&
                mp->atom_array[nbp[i].atoms[i2]].color == 6) {
              continue;
            }
            sum = 0;
            sum += mp->atom_array[nbp[i].atoms[i1]].color *
                   mp->bond_array[nbp[i].bonds[i1]].color;
            sum += mp->atom_array[nbp[i].atoms[i2]].color *
                   mp->bond_array[nbp[i].bonds[i2]].color;
            prod = 1;
            prod *= mp->atom_array[nbp[i].atoms[i1]].color;
            prod &= 0xFFF;
            prod *= mp->atom_array[nbp[i].atoms[i2]].color;
            prod &= 0xFFF;
            seed = NEXT_SEED(seed, sum);
            seed = NEXT_SEED(seed, prod);
            ADD_BIT(fp_counts, ncounts, seed);
            // seed = NEXT_SEED(seed, (sum*prod)&0xFFF);
            // ADD_BIT(fp_counts, ncounts, seed);
            result += 1;
            nmulti = 0;
            if (mp->bond_array[nbp[i].bonds[i1]].color >= 2) {
              nmulti++;
            }
            if (mp->bond_array[nbp[i].bonds[i2]].color >= 2) {
              nmulti++;
            }
            // add bit for hetero-substituted unsaturated hetero atom if
            // degree>1 is well-defined => catch e.g. nitroso vs. nitro
            if (ap->color != 6 && nmulti > 0 &&
                (!as_query || ap->sub_desc == SUB_AS_IS ||
                 (ap->sub_desc != NONE && ap->sub_desc != SUB_MORE &&
                  ap->sub_desc == degree[i] + SUB_ONE - 1))) {
              ADD_BIT(fp_counts, ncounts,
                      NEXT_SEED(seed, nmulti + 173 * degree[i]));
              result += 1;
            }
          }
        }
      }

      if (degree[i] <= 2) {
        continue;
      }
      for (int i1 = 0; i1 < nbp[i].n_ligands; i1++) {
        if (nbp[i].atoms[i1] + 1 == exclude_atom) {
          continue;
        }
        for (int i2 = i1 + 1; i2 < nbp[i].n_ligands; i2++) {
          if (nbp[i].atoms[i2] + 1 == exclude_atom) {
            continue;
          }
          for (int i3 = i2 + 1; i3 < nbp[i].n_ligands; i3++) {
            if (nbp[i].atoms[i3] + 1 == exclude_atom) {
              continue;
            }
            seed = AUGMENTED_ATOM_SEED;
            seed = NEXT_SEED(seed, ap->color);
            if (mp->atom_array[nbp[i].atoms[i1]].color == 0) {
              continue;
            }
            if (mp->atom_array[nbp[i].atoms[i2]].color == 0) {
              continue;
            }
            if (mp->atom_array[nbp[i].atoms[i3]].color == 0) {
              continue;
            }
            if (mp->bond_array[nbp[i].bonds[i1]].color == 0) {
              continue;
            }
            if (mp->bond_array[nbp[i].bonds[i2]].color == 0) {
              continue;
            }
            if (mp->bond_array[nbp[i].bonds[i3]].color == 0) {
              continue;
            }
            /* count hetero neighbours */
            nqtmp = 0;
            if (mp->atom_array[nbp[i].atoms[i1]].color != 6) {
              nqtmp++;
            }
            if (mp->atom_array[nbp[i].atoms[i2]].color != 6) {
              nqtmp++;
            }
            if (mp->atom_array[nbp[i].atoms[i3]].color != 6) {
              nqtmp++;
            }

            /* make sure to add some bits for really odd ones */
            nmulti = 0;
            if (mp->bond_array[nbp[i].bonds[i1]].color == 2) {
              nmulti++;
            }
            if (mp->bond_array[nbp[i].bonds[i1]].color == 3) {
              nmulti++;
            }
            if (mp->bond_array[nbp[i].bonds[i2]].color == 2) {
              nmulti++;
            }
            if (mp->bond_array[nbp[i].bonds[i2]].color == 3) {
              nmulti++;
            }
            if (mp->bond_array[nbp[i].bonds[i3]].color == 2) {
              nmulti++;
            }
            if (mp->bond_array[nbp[i].bonds[i3]].color == 3) {
              nmulti++;
            }

            sum = 0;
            sum += mp->atom_array[nbp[i].atoms[i1]].color *
                   mp->bond_array[nbp[i].bonds[i1]].color;
            sum += mp->atom_array[nbp[i].atoms[i2]].color *
                   mp->bond_array[nbp[i].bonds[i2]].color;
            sum += mp->atom_array[nbp[i].atoms[i3]].color *
                   mp->bond_array[nbp[i].bonds[i3]].color;
            prod = 1;
            prod *= mp->atom_array[nbp[i].atoms[i1]].color;
            prod &= 0xFFF;
            prod *= mp->atom_array[nbp[i].atoms[i2]].color;
            prod &= 0xFFF;
            prod *= mp->atom_array[nbp[i].atoms[i3]].color;
            prod &= 0xFFF;
            seed = NEXT_SEED(seed, sum);
            seed = NEXT_SEED(seed, prod);
            if ((nqtmp > 2 || nmulti >= 2) && ap->color == 6) {
              ADD_BIT(fp_counts, ncounts, seed);
              ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 73));
              result += 2;
            }
            seed = NEXT_SEED(seed, (sum * prod) & 0xFFF);
            ADD_BIT(fp_counts, ncounts, seed);
            result++;
            if (nmulti >= 2 || ap->color > 6)  // make sure R-NO2 is covered
            {
              ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 53));
              result++;
            }
          }
        }
      }
    }
  }

  if (which_bits & USE_AUGMENTED_BOND) {
    /* Set bits for all bonds with both end-degrees > 2 */
    bp = mp->bond_array;
    for (int i = 0; i < mp->n_bonds; i++, bp++) {
      if (bp->atoms[0] == exclude_atom) {
        continue;
      }
      if (bp->atoms[1] == exclude_atom) {
        continue;
      }
      ai1 = bp->atoms[0] - 1;
      ai2 = bp->atoms[1] - 1;
      if (degree[ai1] <= 2) {
        continue;
      }
      if (degree[ai2] <= 2) {
        continue;
      }
      for (int i1 = 0; i1 < nbp[ai1].n_ligands; i1++) {
        for (int i2 = i1 + 1; i2 < nbp[ai1].n_ligands; i2++) {
          /* don't reuse current bond */
          if (nbp[ai1].bonds[i1] == i) {
            continue;
          }
          if (nbp[ai1].bonds[i2] == i) {
            continue;
          }
          if (nbp[ai1].atoms[i1] + 1 == exclude_atom) {
            continue;
          }
          if (nbp[ai1].atoms[i2] + 1 == exclude_atom) {
            continue;
          }
          /* don't use masked atoms/bonds */
          if (mp->atom_array[nbp[ai1].atoms[i1]].color == 0) {
            continue;
          }
          if (mp->atom_array[nbp[ai1].atoms[i2]].color == 0) {
            continue;
          }
          if (mp->bond_array[nbp[ai1].bonds[i1]].color == 0) {
            continue;
          }
          if (mp->bond_array[nbp[ai1].bonds[i2]].color == 0) {
            continue;
          }
          sumi = 0;
          sumi += mp->atom_array[nbp[ai1].atoms[i1]].color *
                  mp->bond_array[nbp[ai1].bonds[i1]].color;
          sumi += mp->atom_array[nbp[ai1].atoms[i2]].color *
                  mp->bond_array[nbp[ai1].bonds[i2]].color;
          prodi = 1;
          prodi *= mp->atom_array[nbp[ai1].atoms[i1]].color +
                   mp->bond_array[nbp[ai1].bonds[i1]].color;
          prodi *= mp->atom_array[nbp[ai1].atoms[i2]].color +
                   mp->bond_array[nbp[ai1].bonds[i2]].color;
          sumi &= 0x0FFF;
          prodi &= 0x0FFF;
          for (int j1 = 0; j1 < nbp[ai2].n_ligands; j1++) {
            for (int j2 = j1 + 1; j2 < nbp[ai2].n_ligands; j2++) {
              /* don't reuse current bond */
              if (nbp[ai2].bonds[j1] == i) {
                continue;
              }
              if (nbp[ai2].bonds[j2] == i) {
                continue;
              }
              if (nbp[ai2].atoms[j1] + 1 == exclude_atom) {
                continue;
              }
              if (nbp[ai2].atoms[j2] + 1 == exclude_atom) {
                continue;
              }
              /* don't use masked atoms/bonds */
              if (mp->atom_array[nbp[ai2].atoms[j1]].color == 0) {
                continue;
              }
              if (mp->atom_array[nbp[ai2].atoms[j2]].color == 0) {
                continue;
              }
              if (mp->bond_array[nbp[ai2].bonds[j1]].color == 0) {
                continue;
              }
              if (mp->bond_array[nbp[ai2].bonds[j2]].color == 0) {
                continue;
              }
              sumj = 0;
              sumj += mp->atom_array[nbp[ai2].atoms[j1]].color *
                      mp->bond_array[nbp[ai2].bonds[j1]].color;
              sumj += mp->atom_array[nbp[ai2].atoms[j2]].color *
                      mp->bond_array[nbp[ai2].bonds[j2]].color;
              prodj = 1;
              prodj *= mp->atom_array[nbp[ai2].atoms[j1]].color +
                       mp->bond_array[nbp[ai2].bonds[j1]].color;
              prodj *= mp->atom_array[nbp[ai2].atoms[j2]].color +
                       mp->bond_array[nbp[ai2].bonds[j2]].color;
              sumj &= 0x0FFF;
              prodj &= 0x0FFF;
              seed = AUGMENTED_BOND_SEED;
              seed = NEXT_SEED(seed, bp->color);
              seed = NEXT_SEED(seed, prodi + prodj);
              seed = NEXT_SEED(seed, sumi * sumj);
              ADD_BIT(fp_counts, ncounts, seed);
              result++;
            }
          }
        }
      }
    }
  }

  if (1 * which_bits & USE_HCOUNT_PAIR) {
    /* generate bits for hydrogen counted described bonds */
    bp = mp->bond_array;
    for (int i = 0; i < mp->n_bonds; i++, bp++) {
      if (H_count[bp->atoms[0]] == 0 && H_count[bp->atoms[1]] == 0) {
        continue;
      }
      if (bp->atoms[0] == exclude_atom) {
        continue;
      }
      if (bp->atoms[1] == exclude_atom) {
        continue;
      }
      /* Don't consider CC single bonds */
      if (mp->atom_array[bp->atoms[0] - 1].color == 6 &&
          mp->atom_array[bp->atoms[1] - 1].color == 6 && bp->color == 1) {
        continue;
      }
      /* Don't consider explicit AH bonds */
      if (mp->atom_array[bp->atoms[0] - 1].color == 0) {
        continue;
      }
      if (mp->atom_array[bp->atoms[1] - 1].color == 0) {
        continue;
      }
      for (int j1 = 0; j1 <= H_count[bp->atoms[0]]; j1++) {
        for (int j2 = 0; j2 <= H_count[bp->atoms[1]]; j2++) {
          if (j1 + j2 == 0) {
            continue; /* at least one hydrogen */
          }
          // unsaturation triggers bit like a hydrogen
          if (!unsaturated[bp->atoms[0] - 1] &&
              !unsaturated[bp->atoms[1] - 1] && j1 * j2 == 0 &&
              bp->color <= 1) {
            continue;
          }
          if (j1 + j2 > 3) {
            continue; /* at most 3 hydrogens */
          }
          // seed = HCOUNT_PAIR_SEED + 53*(j1+j2) + 7*(j1+1)*(j2+1);
          seed = NEXT_SEED(HCOUNT_PAIR_SEED, 7 * (j1 + 1) * (j2 + 1));
          seed = NEXT_SEED(seed, 53 * (j1 + j2));
          seed = NEXT_SEED(seed, bp->color);
          seed = NEXT_SEED(seed, mp->atom_array[bp->atoms[0] - 1].color +
                                     mp->atom_array[bp->atoms[1] - 1].color);
          seed = NEXT_SEED(seed, mp->atom_array[bp->atoms[0] - 1].color *
                                     mp->atom_array[bp->atoms[1] - 1].color);
          ADD_BIT(fp_counts, ncounts, seed);
          result++;
          if (bp->color > 1) {
            seed = NEXT_SEED(seed, 83);
            ADD_BIT(fp_counts, ncounts, seed);
            result++;
            /* Add more bits for really special pairs */
            if (j1 + j2 == 1 &&
                (mp->atom_array[bp->atoms[0] - 1].color != 6 ||
                 mp->atom_array[bp->atoms[1] - 1].color != 6) &&
                bp->color <= 3) {
              seed = NEXT_SEED(seed, 91);
              ADD_BIT(fp_counts, ncounts, seed);
              seed = NEXT_SEED(seed, 97);
              ADD_BIT(fp_counts, ncounts, seed);
              seed = NEXT_SEED(seed, 103);
              ADD_BIT(fp_counts, ncounts, seed);
              result += 3;
            }
          }
        }
      }
    }
  }

  if (which_bits & USE_HCOUNT_PATH) {
    /* generate a short path for each atom that has a hydrogen */
    seed = HCOUNT_PATH_SEED;
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (H_count[i + 1] == 0) {
        continue;
      }
      /* don't consider non methyl carbon atoms */
      if (ap->color == 6 && H_count[i + 1] < 3 && degree[i] < 3) {
        continue;
      }

      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      if (ap->color != 6) {
        result +=
            SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices,
                           1, 2, 5, /* path length 2 to 4 */  // EVG
                           i, 0, -1,
                           IGNORE_PATH_SYMBOL |
                               // IGNORE_TERM_SYMBOL |
                               PROCESS_CHAINS,
                           exclude_atom);
      } else {
        if (degree[i] > 2)  // tertiary hydrogen
        {
          if (0) {  // Class disabled to save bit density
            result += SetPathBitsRec(
                mp, nbp, fp_counts, ncounts, NEXT_SEED(seed, 101),
                touched_indices, 1, 2, 2, /* path length 2 to 3 */
                i, 0, -1,
                IGNORE_PATH_SYMBOL | IGNORE_TERM_SYMBOL | PROCESS_CHAINS,
                exclude_atom);
          }
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, NEXT_SEED(seed, 103),
              touched_indices, 1, 2, 5, /* path length 2 to 3 */
              i, 0, -1,
              IGNORE_PATH_SYMBOL | FORCED_HETERO_END | /* i-Pr...Q */
                  PROCESS_CHAINS,
              exclude_atom);
        }
        if (H_count[i + 1] >= 3)  // methyl
        {
          if (0) {  // Class disabled to save bit density
            result += SetPathBitsRec(
                mp, nbp, fp_counts, ncounts, NEXT_SEED(seed, 1103),
                touched_indices, 1, 2, 3, /* path length 2 to 3 */
                i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
          }
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, NEXT_SEED(seed, 1105),
              touched_indices, 1, 3, 6, /* path length 4 to 4 */
              i, 0, -1,
              IGNORE_PATH_SYMBOL | FORCED_HETERO_END | /* Me...Q */
                  PROCESS_CHAINS,
              exclude_atom);
        }
      }
      touched_indices[i] = 0; /* down-dating */
      seed = old_seed;
      if (ap->color == 6) {
        continue;
      }

      if (H_count[i + 1] > 1) /* catch the difference between NH and NH2! */
      {
        touched_indices[i] = 1; /* updating */
        old_seed = seed;
        seed = NEXT_SEED(seed, 113);
        seed = NEXT_SEED(seed, ap->color);
        result +=
            SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices,
                           1, 1, 5, /* path length 1 to 5 */
                           i, 0, -1,
                           // FORCED_HETERO_END |
                           IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
        touched_indices[i] = 0; /* down-dating */
        seed = old_seed;
      }

      for (int j = 1; j < H_count[i + 1]; j++) {
        seed = NEXT_SEED(seed, 61 * j);
        seed = NEXT_SEED(seed, ap->color);
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
      }
      seed = old_seed;
    }
  }

  /* Compute ring paths */
  ap = mp->atom_array;
  for (int i = 0; i < mp->n_atoms; i++, ap++) {
    if (atom_status[i] <= 0) {
      ap->color = 0;
    }
    if (i + 1 == exclude_atom) {
      ap->color = 0;
    }
    if (ap->color == 0) {
      continue;
    }
  }

  /* remove all bonds with only non-ring atoms from consideration */
  bp = mp->bond_array;
  for (int i = 0; i < mp->n_bonds; i++, bp++) {
    if (bond_status[i] <= 0 && atom_status[bp->atoms[0] - 1] == 0 &&
        atom_status[bp->atoms[1] - 1] == 0) {
      bp->color = 0;
    }
    if (bp->atoms[0] == exclude_atom) {
      bp->color = 0;
    }
    if (bp->atoms[1] == exclude_atom) {
      bp->color = 0;
    }
  }

  if (which_bits & USE_RING_PATH) {
    seed = RING_PATH_SEED;
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      if (0) {  // Class disabled to save bit density
        result +=
            SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices,
                           1, 2, 3, i, 0, -1, PROCESS_CHAINS, exclude_atom);
      }

      if (ap->color > 5 && ap->color < 10 &&
          atom_status[i] > 2)  // only start at common light atoms
      {
        seed = NEXT_SEED(seed, 61);
        result +=
            SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices,
                           1, 3, 8, /* 3 to 8 */
                           i, 0, -1,
                           STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                               PROCESS_RING_CLOSURES,
                           exclude_atom);
      }
      seed = old_seed;

      touched_indices[i] = 0; /* down-dating */
    }
  }

  /* Set the color property to represent all different atom types */
  ap = mp->atom_array;
  for (int i = 0; i < mp->n_atoms; i++, ap++) {
    ap->color = AtomicNumberFromSymbol(ap->atom_symbol);
    if (ap->color == 1 || i + 1 == exclude_atom) {
      ap->color = 0; /* ignore hydrogens */
    } else {
      ap->color = ANY_COLOR; /* treat all other atoms alike */
    }
  }

  /* Set the color property to represent the different bond type classes */
  bp = mp->bond_array;
  for (int i = 0; i < mp->n_bonds; i++, bp++) {
    if (bp->bond_type == SINGLE) {
      bp->color = 1;
    } else if (bp->bond_type == DOUBLE) {
      bp->color = 2;
    } else if (bp->bond_type == TRIPLE) {
      bp->color = 3;
    } else if (bp->bond_type == AROMATIC) {
      bp->color = 4;
    } else {
      bp->color = 0;
    }
    if (bp->atoms[0] == exclude_atom) {
      bp->color = 0;
    }
    if (bp->atoms[1] == exclude_atom) {
      bp->color = 0;
    }
  }

  // Here, we have all non-trivial atoms mapped to ANY_COLOR while the
  // bond type is retained.
  if (which_bits & USE_BOND_PATH) {
    seed = BOND_PATH_SEED;
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      // start at branch node on ring
      if (degree[i] > 2 && atom_status[i] > 1) {
        result += SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed,
                                 touched_indices, 1, 4, 4, /* was 4 to 4 */
                                 i, 0, -1, PROCESS_CHAINS, exclude_atom);
        if (ap->rsize_flags & SPECIAL_RING) {
          result +=
              SetPathBitsRec(mp, nbp, fp_counts, ncounts, NEXT_SEED(seed, 217),
                             touched_indices, 1, 5, 5, i, 0, -1,
                             IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
        }
      }
      /* Add other bits to catch poorly specified ring closures */
      seed = old_seed;
      seed = NEXT_SEED(seed, 11);
      result += SetPathBitsRec(
          mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 4,
          6, /* was 5 to 6 */
          i, 0, -1,
          // DEBUG_PATH |
          FORCED_RING_PATH | IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
          exclude_atom);
      /* Add bits for paths starting with rare bond orders */
      for (int j = 0; j < nbp[i].n_ligands; j++) {
        bp = &mp->bond_array[nbp[i].bonds[j]];
        ai = nbp[i].atoms[j];
        if (ai + 1 == exclude_atom) {
          continue;
        }
        if (bp->color == 0) {
          continue;
        }
        if (bp->bond_type != DOUBLE && bp->bond_type != TRIPLE) {
          continue;
        }
        if (atom_status[i] <= 0 && bp->bond_type != TRIPLE) {
          continue;
        }
        seed = old_seed;
        seed = NEXT_SEED(seed, bp->color * 413);
        touched_indices[ai] = 1; /* updating */
        result += SetPathBitsRec(
            mp, nbp, fp_counts, ncounts, seed, touched_indices, 2, 4, 5, ai, 0,
            i, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES | PROCESS_CHAINS,
            exclude_atom);
        touched_indices[ai] = 0; /* down-dating */
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
  }

  /* Set the color property to represent the different atom type classes */
  ap = mp->atom_array;
  for (int i = 0; i < mp->n_atoms; i++, ap++) {
    if (i + 1 == exclude_atom) {
      ap->color = 0;
      continue;
    }
    {
      if (ap->atom_symbol == "H") {
        ap->color = 0; /* ignore hydrogens */
      } else if (ap->atom_symbol == "D") {
        ap->color = 0; /* ignore hydrogens */
      } else if (ap->atom_symbol == "T") {
        ap->color = 0; /* ignore hydrogens */
      } else if (ap->atom_symbol == "C") {
        ap->color = 6; /* carbon second row elements are one class */
      } else if (ap->atom_symbol == "N") {
        ap->color = 8; /* nitrogen, oxigen, and sulfur are one class */
      } else if (ap->atom_symbol == "O") {
        ap->color = 8; /* nitrogen, oxigen, and sulfur are one class */
      } else if (ap->atom_symbol == "S") {
        ap->color = 8; /* nitrogen, oxigen, and sulfur are one class */
      } else if (ap->atom_symbol == "Q") {
        ap->color = 8; /* nitrogen, oxigen, and sulfur are one class */
      } else if (ap->atom_symbol == "A") {
        ap->color = 0;
      } else {
        tmp = AtomicNumberFromSymbol(ap->atom_symbol);
        if (1 < tmp && tmp < 115) {
          ap->color = 8; /* non-carbon is the only second class */
        } else {         /* This could be R atoms or other odd things */
          ap->color = 0;
        }
      }
    }
    if (ap->color > 115) {
      ap->color = 0; /* ignore special atom types */
    }
  }

  // Here, we unify atom types to carbon/hetero distinction
  if (which_bits & USE_HCOUNT_CLASS_PATH) {
    /* generate a short path for each atom that has a hydrogen */
    seed = HCOUNT_CLASS_PATH_SEED;
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (H_count[i + 1] == 0) {
        continue;
      }
      if (ap->color == 6 && H_count[i + 1] < 2) {
        continue;
      }
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      touched_indices[i] = 1; /* updating */
      if (ap->color == 6) {
        result += SetPathBitsRec(
            mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 2,
            4, /* path length 1 to 3 */
            i, 0, -1, FORCED_HETERO_END | PROCESS_CHAINS, exclude_atom);
      } else {
        result += SetPathBitsRec(
            mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 2, 5, i, 0,
            -1, IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS,
            exclude_atom);
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
  }

  /* Do a small first pass with defined bond types */
  if (which_bits & USE_ATOM_CLASS_PATH) {
    // seed = ATOM_CLASS_PATH_SEED+117;
    seed = NEXT_SEED(ATOM_CLASS_PATH_SEED, 117);
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (ap->color == 6 && degree[i] < 3) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      if (ap->color == 6) {
        if (0) {  // Class disabled to save bit density
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
              4, /* path length 3 to 4 */
              i, 0, -1, FORCED_HETERO_END | IGNORE_PATH_SYMBOL | PROCESS_CHAINS,
              exclude_atom);
        }
      } else {
        if (0) {  // Class disabled to save bit density
          result +=
              SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices,
                             1, 2, 3, /* path length 3 to 3 */
                             i, 0, -1,
                             // IGNORE_PATH_SYMBOL |
                             PROCESS_CHAINS, exclude_atom);
        }
        if (0) {  // Class disabled to save bit density
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 4,
              5, /* path length 4 to 5 */
              i, 0, -1, IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS,
              exclude_atom);
        }
        result += SetPathBitsRec(
            mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
            4, /* path length 4 to 7 */
            i, 0, -1, FORCED_RING_PATH | PROCESS_RING_CLOSURES | PROCESS_CHAINS,
            exclude_atom);
      }
      touched_indices[i] = 0;
      seed = old_seed;
    }
  }

  /* Set the color property to only a single class */
  bp = mp->bond_array;
  for (int i = 0; i < mp->n_bonds; i++, bp++) {
    if (SINGLE <= bp->bond_type && bp->bond_type <= ANY_BOND) {
      bp->color = 5;
    } else {
      bp->color = 0;
    }
    if (bp->atoms[0] == exclude_atom) {
      bp->color = 0;
    }
    if (bp->atoms[1] == exclude_atom) {
      bp->color = 0;
    }
  }

  // Here, we've unified atom types to carbon/hetero and made all bond types
  // identical
  if (which_bits & USE_ATOM_CLASS_PATH) {
    seed = ATOM_CLASS_PATH_SEED;
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (ap->color == 6) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      // if (ap->color == 6)
      {
        if (0 * degree[i] > 2) {  // Class disabled to save bit density
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
              3, /* path length 3 to 3 */
              i, 0, -1, FORCED_HETERO_END | PROCESS_CHAINS, exclude_atom);
        }
      }
      // else
      {
        if (0) {  // Class disabled to save bit density
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
              4, /* path length 3 to 4 */
              i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
        }
        seed = NEXT_SEED(seed, 23 + ap->color * 19);
        result += SetPathBitsRec(
            mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
            9, /* path length 3 to 9 */
            i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES, exclude_atom);
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
    /* Q-Q and Q-C ring bond count */
    qq_count = 0;
    qc_count = 0;
    bp = mp->bond_array;
    for (int i = 0; i < mp->n_bonds; i++, bp++) {
      if (bond_status[i] == 0) {
        continue;
      }
      if (bp->color == 0) {
        continue;
      }
      if (bp->atoms[0] == exclude_atom) {
        continue;
      }
      if (bp->atoms[1] == exclude_atom) {
        continue;
      }
      ai1 = mp->atom_array[bp->atoms[0] - 1].color;
      ai2 = mp->atom_array[bp->atoms[1] - 1].color;
      if (ai1 == 0 || ai2 == 0) {
        continue;
      }
      if (ai1 == 6 && ai2 == 6) {
        continue;
      }
      if (ai1 != 6 && ai2 != 6) {
        qq_count++;
        for (int j = 3; j < 9; j++) { /* set bits for not too large ring size */
          if (bp->rsize_flags & (1 << j)) {
            ADD_BIT(fp_counts, ncounts,
                    NEXT_SEED(ATOM_CLASS_PATH_SEED * 17, j * 8));
            result++;
          }
        }
      } else {
        qc_count++;
        for (int j = 3; j < 9; j++) { /* set bits for not too large ring size */
          if (bp->rsize_flags & (1 << j)) {
            ADD_BIT(fp_counts, ncounts,
                    NEXT_SEED(ATOM_CLASS_PATH_SEED * 19, j * 8));
            result++;
          }
        }
      }
    }
    seed = 2 * ATOM_CLASS_PATH_SEED + 3;
    for (int i = 1; i <= qq_count; i = (int)(1 + i * 1.5)) {
      seed = NEXT_SEED(seed, i * 153);
      ADD_BIT(fp_counts, ncounts, seed);
      result++;
      seed = NEXT_SEED(seed, 53);
      if (i <= 1) {
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
      }
    }
    seed = 3 * ATOM_CLASS_PATH_SEED + 5;
    for (int i = 1; i <= MIN(qc_count, 2); i++) {
      seed = NEXT_SEED(seed, i * 157);
      ADD_BIT(fp_counts, ncounts, seed);
      result++;
    }
    for (int i = 3; i <= qc_count; i = (int)(i * 1.8)) {
      seed = NEXT_SEED(seed, i * 157);
      ADD_BIT(fp_counts, ncounts, seed);
      result++;
    }
  }

  /* Compute ring patters with at least one cycle */
  /* remove bonds from consideration that don't have at least one ring atom */
  bp = mp->bond_array;
  for (int i = 0; i < mp->n_bonds; i++, bp++) {
    if (bp->atoms[0] == exclude_atom) {
      continue;
    }
    if (bp->atoms[1] == exclude_atom) {
      continue;
    }
    bp->color = 5;
    /* ignore non-ring bonds */
    if (atom_status[bp->atoms[0] - 1] == 0 &&
        atom_status[bp->atoms[1] - 1] == 0) {
      bp->color = 0;
    } else {
      // may be redundant since atom types have already been unified
      ap = &mp->atom_array[bp->atoms[0] - 1];
      if (ap->color != 0) {
        if (ap->color != 6) {
          ap->color = 8; /* non-carbons in one class */
        }
        if (ap->rsize_flags == 0) {
          ap->color = 0;
        }
      }
      // may be redundant since atom types have already been unified
      ap = &mp->atom_array[bp->atoms[1] - 1];
      if (ap->color != 0) {
        if (ap->color != 6) {
          ap->color = 8; /* non-carbons in one class */
        }
        if (ap->rsize_flags == 0) {
          ap->color = 0;
        }
      }
    }
  }

  // Here, we have bond type ignored and atom types mapped to C and Q
  if (which_bits & USE_RING_PATTERN) {
    /* first process ring bond paths with atom classes */
    seed = RING_PATTERN_SEED;
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      /* Don't process fragments starting at carbon in 6-ring only */
      if (ap->color == 6 && 0 == (ap->rsize_flags & SPECIAL_RING)) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      result +=
          SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices, 1,
                         3, 3, /* ring bond path size 3 to 3 */
                         i, 0, -1, PROCESS_CHAINS, exclude_atom);
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }

    /* Now, we only include complete rings but ignore atom-type */
    /* 'A' atoms are now included nodes */
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      /* add 'A' atom to standard class */
      if (ap->atom_symbol == "A") {
        ap->color = 9;
      }
      if (i + 1 == exclude_atom) {
        ap->color = 0;
      }
      if (ap->color == 0) {
        continue;
      }
      ap->color = 9; /* all ring atoms in same class */
    }
    bp = mp->bond_array;
    for (int i = 0; i < mp->n_bonds; i++, bp++) {
      if (bp->atoms[0] == exclude_atom) {
        bp->color = 0;
      }
      if (bp->atoms[1] == exclude_atom) {
        bp->color = 0;
      }
      if (bond_status[i] <= 0) {
        bp->color = 0;
      }
    }
    // seed = RING_PATTERN_SEED+23;
    seed = NEXT_SEED(RING_PATTERN_SEED, 23);
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color == 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      if (0) {  // Class disabled to save bit density
        result += SetPathBitsRec(
            mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 4,
            17, /* ring size 4 to 17 */
            i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES, exclude_atom);
      }

      seed = old_seed;
      if (0) {                   // Class disabled to save bit density
        if (atom_status[i] > 2)  // start at ring fusion
        {
          seed = NEXT_SEED(seed, 61);
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 6,
              17, /* ring path size 6 to 17 */
              i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
              exclude_atom);
        }
      }

      seed = old_seed;
      if (0) {  // Class disabled to save bit density
        if (degree[i] > 2 && atom_status[i] >= 2)  // start at ring substituents
        {
          seed = NEXT_SEED(seed, 67);
          result += SetPathBitsRec(
              mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 6,
              17, /* ring path size 6 to 17 */
              i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
              exclude_atom);
        }
      }
      seed = old_seed;

      touched_indices[i] = 0; /* down-dating */
    }
  }

  if (which_bits & USE_RING_SIZE_COUNTS) {
    for (int j = 3; j < 10; j++) /* loop through ring_sizes */
    {
      nrbonds = 0;
      bp = mp->bond_array;
      for (int i = 0; i < mp->n_bonds; i++, bp++) {
        if (bp->atoms[0] != exclude_atom && bp->atoms[1] != exclude_atom &&
            (bp->rsize_flags & (1 << j))) {
          nrbonds++;
        }
      }
      seed = RING_SIZE_SEED;
      seed = NEXT_SEED(seed, j * 13);
      for (int i = 1; i < 100; i *= 2) {
        if (nrbonds >= j * i) {
          seed = NEXT_SEED(seed, i);
          /* don't set a bit if just one 5- or 6-ring */
          if ((j != 6 && j != 5) || nrbonds > j * i) {
            ADD_BIT_COUNT(fp_counts, ncounts, seed, nrbonds);
            result++;
          }
          if (j != 6 && j != 5) /* set one more bit for odd rings */
          {
            seed = NEXT_SEED(seed, 17);
            ADD_BIT_COUNT(fp_counts, ncounts, seed, nrbonds);
            result++;
          }
        } else {
          break;
        }
      }
    }
    /* Set bits for different ring sizes connected by a bond */
    for (int j = 0; j < 15; j++) {   /* loop through ring_sizes */
      for (int k = 0; k < 15; k++) { /* loop through ring_sizes */
        rscounts[j][k] = 0;
      }
    }
    /* count connections */
    bp = mp->bond_array;
    for (int i = 0; i < mp->n_bonds; i++, bp++) {
      if (bp->atoms[0] == exclude_atom) {
        continue;
      }
      if (bp->atoms[1] == exclude_atom) {
        continue;
      }
      for (int j = 3; j < 15; j++) { /* loop through ring_sizes */
        for (int k = 3; k < 15; k++) /* loop through ring_sizes */
        {
          if (j == k) {
            continue;
          }
          if ((mp->atom_array[bp->atoms[0] - 1].rsize_flags & (1 << j)) &&
              (mp->atom_array[bp->atoms[1] - 1].rsize_flags & (1 << k))) {
            rscounts[j][k]++;
            rscounts[k][j]++;
          }
        }
      }
    }
    /* set bits */
    for (int j = 3; j < 9; j++) {     /* loop through not too large ring_sizes */
      for (int k = j + 1; k < 9; k++) /* loop through not too large ring_sizes */
      {
        if (rscounts[j][k] == 0) {
          continue;
        }
        seed = 2 * RING_SIZE_SEED;
        seed = NEXT_SEED(seed, (j + k) * 11);
        seed = NEXT_SEED(seed, j * k * 13);
        seed = NEXT_SEED(seed, 19);
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
        /* more special links => more bits */
        if (rscounts[j][k] == 1 || j == 6 || k == 6) {
          continue;
        }
        seed = 2 * RING_SIZE_SEED + 23;
        seed = NEXT_SEED(seed, (j + k) * 11);
        seed = NEXT_SEED(seed, j * k * 13);
        seed = NEXT_SEED(seed, 19);
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
      }
    }
  }

  /* Set the color property to represent all different atom types */
  ap = mp->atom_array;
  for (int i = 0; i < mp->n_atoms; i++, ap++) {
    ap->color = AtomicNumberFromSymbol(ap->atom_symbol);
    if (ap->color <= 1) {
      ap->color = 0; /* ignore hydrogens */
    }
    /* mark special atom types */
    if (ap->color > 115) {
      ap->color = -1;
    }
    if (ap->atom_symbol == "A") {
      ap->color = -1;
    }
    if (i + 1 == exclude_atom) {
      ap->color = 0;
    }
    if (ap->color > 1 && (!as_query || ap->sub_desc == SUB_AS_IS ||
                          (ap->sub_desc != NONE && ap->sub_desc != SUB_MORE &&
                           ap->sub_desc == degree[i] + SUB_ONE - 1))) {
      ap->color += 32 * degree[i];
    } else {
      ap->color = 0;
    }
  }
  bp = mp->bond_array;
  for (int i = 0; i < mp->n_bonds; i++, bp++) {
    if (SINGLE <= bp->bond_type && bp->bond_type <= ANY_BOND) {
      bp->color = 5;
    } else {
      bp->color = 0;
    }
    if (bp->atoms[0] == exclude_atom) {
      bp->color = 0;
    }
    if (bp->atoms[1] == exclude_atom) {
      bp->color = 0;
    }
  }

  /*
   * add special bits for paths starting at atoms with >= 4 neighbours or
   * methyl atoms
   *
   * atom color represent degree bond color identical for all bonds to non-H
   */
  if (which_bits & USE_DEGREE_PATH) {
    seed = DEGREE_PATH_SEED;

    // set bits for degree paths starting with special carbon atoms
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (ap->color <= 0) {
        continue;  // only process if degree defined
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      /* don't start on usual atoms */
      if (degree[i] <= 3 &&       // include if high degree
          atom_status[i] <= 2 &&  // include if ring fusion
          degree[i] != 1) {       // include if terminal atom
        continue;
      }
      /* only keep terminals if methyl */
      if (degree[i] == 1 && ap->atom_symbol != "C") {
        continue;
      }
      // ap->atom_symbol, i+1, degree[i], atom_status[i]);
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, ap->color);
      // ap->atom_symbol, i+1, degree[i], ap->sub_desc);
      result +=
          SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices, 1,
                         2, 4, /* path length 1 to 3 */
                         i, 0, -1,
                         // IGNORE_TERM_SYMBOL |
                         IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
      /* special CH fusion atoms */
      if (atom_status[i] > 2 && H_count[i + 1] >= 1) {
        seed = NEXT_SEED(seed, 219);
        result += SetPathBitsRec(
            mp, nbp, fp_counts, ncounts, seed, touched_indices, 1, 2, 5, i, 0,
            -1, IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }

    // set bits for degree paths starting with hetero atoms
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      tmp = AtomicNumberFromSymbol(ap->atom_symbol);
      if (1 >= tmp || tmp >= 115) {
        continue;
      }
      if (tmp == 6) {
        continue;
      }
      if (tmp < 10) {
        continue;  // exclude common hetero atoms
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, tmp);
      // ap->atom_symbol, i+1, degree[i], ap->sub_desc);
      result +=
          SetPathBitsRec(mp, nbp, fp_counts, ncounts, seed, touched_indices, 1,
                         2, 2, i, 0, -1, PROCESS_CHAINS, exclude_atom);
      if (0) {  // might overly populate complexes
        result += SetPathBitsRec(
            mp, nbp, fp_counts, ncounts, 101 + seed, touched_indices, 1, 2, 2,
            i, 0, -1, FORCED_RING_PATH | PROCESS_CHAINS, exclude_atom);
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
  }

  if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS)) {
    /* Collect length_matrix */
    /* allocate storage length_matrix */
    length_matrix.assign(mp->n_atoms, std::vector<int>(mp->n_atoms, 0));
    for (int i = 0; i < mp->n_atoms; i++) {
      touched_indices[i] = 0;
    }
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      // if (FALSE) fprintf(stderr, "starting path search at atom %d(%d)\n",
      // i+1, ap->color);
      SetPathLengthFlags(mp, touched_indices, i, 0, i,
                         12, /* path perception distance <= 12 */
                         length_matrix, nbp, exclude_atom);
      touched_indices[i] = 0; /* down-dating */
    }
  }

  /*
   * This screen class will catch non-linear fragments composed of rather
   * frequent linear sub-fragments
   */
  if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS)) {
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      ap->color = AtomicNumberFromSymbol(ap->atom_symbol);
      if (ap->atom_symbol == "H") {
        ap->color = 0; /* ignore hydrogens */
      } else if (ap->atom_symbol == "D") {
        ap->color = 0; /* ignore hydrogens */
      } else if (ap->atom_symbol == "T") {
        ap->color = 0; /* ignore hydrogens */
      } else if (ap->atom_symbol == "Q") {
        ap->color = HETERO;
      } else if (ap->atom_symbol == "A") {
        ap->color = GENERIC;
      } else if (ap->atom_symbol == "L") {
        ap->color = GENERIC;
      } else if (ap->atom_symbol == "C") {
        ap->color = 6; /* carbon second row elements are one class */
        if (cdegree[i] >= 3) {
          ap->color = CSP3;
        }
      } else if (ap->color > 1 && ap->color < 115) {
        ap->color = HETERO;
      } else { /* This could be R atoms or other odd things */
        ap->color = 0;
      }
      if (i + 1 == exclude_atom) {
        ap->color = 0;
      }
    }
    /*
     * Bond colors are already set OK, i.e. equal for A-H bonds
     */
    /* NOP */

    /* Now we start setting bits */
    ap = mp->atom_array;
    for (int i = 0; i < mp->n_atoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      /* Spiders have at least three legs (;-) */
      if (degree[i] < 3) {
        continue;
      }
      /* Spider needs to be special atom or carbon */
      if (ap->color != CSP3 && ap->color != 6) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      std::fill(csp3.begin(), csp3.end(), 0);
      std::fill(hetero.begin(), hetero.end(), 0);
      if (which_bits & USE_CLASS_SPIDERS) {
        SpecialNeighboursRec(mp, touched_indices, 1, i, MAX_SPIDER, csp3.data(),
                             hetero.data(), nbp, exclude_atom);
      }
      touched_indices[i] = 0; /* down-dating */

      /* set bits for spiders with one CSP3 atom and two heteros */
      if (which_bits & USE_CLASS_SPIDERS) {
        for (int j = 1; j <= MAX_SPIDER; j++) {
          if (csp3[j] == 0) {
            continue;
          }
          seed = CLASS_SPIDER_SEED;
          if (ap->color == HETERO) {
            // seed = CLASS_SPIDER_SEED+HETERO*8+CSP3*11;
            seed = NEXT_SEED(seed, HETERO * 8);
            seed = NEXT_SEED(seed, CSP3 * 11);
          } else {
            // seed = CLASS_SPIDER_SEED+6*8+CSP3*11;
            seed = NEXT_SEED(seed, 6 * 8);
            seed = NEXT_SEED(seed, CSP3 * 11);
          }
          for (int j1 = 1; j1 <= MAX_SPIDER; j1++) {
            tmp1 = hetero[j1];
            if (tmp1 <= 0) {
              continue;
            }
            for (int j2 = j1; j2 <= MAX_SPIDER; j2++) {
              tmp2 = hetero[j2];
              if (j2 == j1) {
                tmp2--; /* consumed in outer loop */
              }
              if (tmp2 <= 0) {
                continue;
              }
              old_seed = seed;
              seed = NEXT_SEED(seed, j);
              seed = NEXT_SEED(seed, j1 + j2);
              ADD_BIT(fp_counts, ncounts, seed);
              result++;
              seed = old_seed;
            }
          }
        }
      }

      /* don't hetero-spider normal carbons */
      if (ap->color != CSP3) {
        continue;
      }
      /* set bits for spiders with three defined HETERO atoms */
      if (which_bits & USE_CLASS_SPIDERS) {
        for (int j = 1; j <= MAX_SPIDER; j++) {
          if (hetero[j] == 0) {
            continue;
          }
          seed = CLASS_SPIDER_SEED;
          if (ap->color == HETERO) {
            // seed = CLASS_SPIDER_SEED+HETERO*8+HETERO*11;
            seed = NEXT_SEED(seed, HETERO * 8);
            seed = NEXT_SEED(seed, HETERO * 11);
          } else {
            // seed = CLASS_SPIDER_SEED+6*8+HETERO*11;
            seed = NEXT_SEED(seed, 6 * 8);
            seed = NEXT_SEED(seed, HETERO * 11);
          }
          for (int j1 = j; j1 <= MAX_SPIDER; j1++) {
            tmp1 = hetero[j1];
            if (j1 == j) {
              tmp1--; /* we've consumed this one in outer loop */
            }
            if (tmp1 <= 0) {
              continue;
            }
            for (int j2 = j1; j2 <= MAX_SPIDER; j2++) {
              tmp2 = hetero[j2];
              if (j2 == j) {
                tmp2--; /* consumed in outer loop */
              }
              if (j2 == j1) {
                tmp2--; /* consumed in outer loop */
              }
              if (tmp2 <= 0) {
                continue;
              }
              old_seed = seed;
              seed = NEXT_SEED(seed, j * j1 * j2);
              ADD_BIT(fp_counts, ncounts, seed);
              result++;
              // Additional bits for quarternary centers
              if (degree[i] > 3) {
                ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 501));
                result++;
              }
              seed = old_seed;
            }
          }
        }
      }
    }

    /**
     * Collect bits that represent feature/path_length/feature triples.
     */
    if (which_bits & USE_FEATURE_PAIRS) {
      /* set feature flags in atom colors */
      ap = mp->atom_array;
      for (int i = 0; i < mp->n_atoms; i++, ap++) {
        if (ap->atom_symbol == "C") {
          flags = C_FLAG;
        } else if (ap->atom_symbol == "O") {
          flags = O_FLAG;
        } else if (ap->atom_symbol == "N") {
          flags = N_FLAG;
        } else if (ap->atom_symbol == "S") {
          flags = S_FLAG;
        } else if (ap->atom_symbol == "P") {
          flags = P_FLAG;
        } else if (AtomSymbolMatch(ap->atom_symbol, "F,Cl,Br,I,At")) {
          flags = X_FLAG;
        } else {
          flags = 0;
        }
        if (ap->color == HETERO) {
          flags |= HETERO_FLAG;
        }
        if (cdegree[i] >= 3) {
          flags |= CSP3_FLAG;
        }
        if (degree[i] >= 4) {
          flags |= QUART_FLAG;
        }
        if (atom_status[i] > 0 && degree[i] >= 3) {
          flags |= RING_SUBST_FLAG;
          if (0 != (ap->rsize_flags & SPECIAL_RING)) {
            flags |= RS_SPECIAL_FLAG;
          }
        }
        if (i + 1 == exclude_atom) {
          flags = 0;
        }
        ap->color = flags;
      }
      /* collect bits for selected feature pairs */
      if (0) {  // Class disabled to save bit density
        result += SetFeatureBits(mp, fp_counts, ncounts, CSP3_FLAG,
                                 HETERO_FLAG, /* to hetero or ring subst */
                                 2, 3,        /* with path length 1 to 9 */
                                 FALSE,       /* don't use count */
                                 TRUE,        /* use atom type flags */
                                 length_matrix, 1237, /* seed = 1237 */
                                 exclude_atom);
      }
      if (0) {  // Class disabled to save bit density
        result += SetFeatureBits(mp, fp_counts, ncounts,
                                 HETERO_FLAG, /* from ring substitution */
                                 HETERO_FLAG, /* to hetero or ring subst */
                                 1, 12,       /* with path length 1 to 10 */
                                 TRUE,        /* use count */
                                 TRUE,        /* use atom type flags */
                                 length_matrix, 1237, /* seed = 1237 */
                                 exclude_atom);
      }
      if (1) {
        result += SetFeatureBits(mp, fp_counts, ncounts,
                                 RING_SUBST_FLAG, /* from ring substitution */
                                 RING_SUBST_FLAG, /* to ring substitution */
                                 5, 7,            /* with path length 1 to 12 */
                                 FALSE,           /* don't use count */
                                 TRUE,            /* use atom type flags */
                                 length_matrix, 2237, /* seed = 2237 */
                                 exclude_atom);
      }
      if (0) {  // Class disabled to save bit density
        result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           RS_SPECIAL_FLAG, /* from special ring substitution */
                           HETERO_FLAG,     /* to hetero atom */
                           2, 4,            /* with path length 1 to 5 */
                           TRUE,            /* don't use count */
                           FALSE,           /* use atom type flags */
                           length_matrix, 3237, /* seed = 3237 */
                           exclude_atom);
      }
      if (1) {
        result += SetFeatureBits(mp, fp_counts, ncounts,
                                 QUART_FLAG,  /* from quartenary atom */
                                 HETERO_FLAG, /* to hetero or ring subst */
                                 1, 8,        /* with path length 1 to 8 */
                                 FALSE,       /* don't use count */
                                 TRUE,        /* use atom type flags */
                                 length_matrix, 4237, /* seed = 4237 */
                                 exclude_atom);
      }
      if (1) {
        result += SetFeatureBits(
            mp, fp_counts, ncounts, QUART_FLAG, /* from quartenary atom */
            RING_SUBST_FLAG, 1, 6,              /* with path length 1 to 8 */
            FALSE,                              /* don't use count */
            TRUE,                               /* use atom type flags */
            length_matrix, 5237,                /* seed = 4237 */
            exclude_atom);
      }
      if (1) {
        result += SetFeatureBits(mp, fp_counts, ncounts,
                                 X_FLAG,          /* from halogen atom */
                                 CSP3_FLAG, 1, 1, /* with path length 1 to 1 */
                                 FALSE,           /* don't use count */
                                 TRUE,            /* use atom type flags */
                                 length_matrix, 15237, /* seed = 15237 */
                                 exclude_atom);
      }
      if (0) {  // Class disabled to save bit density
        result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           RS_SPECIAL_FLAG, /* from special ring substitution */
                           RING_SUBST_FLAG, /* to hetero atom */
                           2, 4,            /* with path length 1 to 5 */
                           TRUE,            /* don't use count */
                           FALSE,           /* use atom type flags */
                           length_matrix, 6237, /* seed = 6237 */
                           exclude_atom);
      }
      if (0) {  // Class disabled to save bit density
        result += SetFeatureBits(mp, fp_counts, ncounts,
                                 HETERO_FLAG,     /* from hetero */
                                 RING_SUBST_FLAG, /* to ring substitution */
                                 1, 6,            /* with path length 1 to 8 */
                                 FALSE,           /* don't use count */
                                 TRUE,            /* use atom type flags */
                                 length_matrix, 7237, /* seed = 4237 */
                                 exclude_atom);
      }

      /* Set bits for ring-subst/ring-subst/hetero triples */
      if (0)  // too many spurious bits
      {
        ap1 = mp->atom_array;
        for (int i1 = 0; i1 < mp->n_atoms; i1++, ap1++) {
          if (i1 + 1 == exclude_atom) {
            continue;
          }
          if (0 == (ap1->color & RING_SUBST_FLAG)) {
            continue;
          }
          /* first atom must be in a non-sixmembered ring */
          if (!(ap1->rsize_flags & SPECIAL_RING)) {
            continue;
          }
          ap2 = mp->atom_array;
          for (int i2 = 0; i2 < mp->n_atoms; i2++, ap2++) {
            if (i1 == i2) {
              continue;
            }
            if (i2 + 1 == exclude_atom) {
              continue;
            }
            if (0 == (ap2->color & RING_SUBST_FLAG)) {
              continue;
            }
            ap3 = mp->atom_array;
            for (int i3 = 0; i3 < mp->n_atoms; i3++, ap3++) {
              if (i1 == i3) {
                continue;
              }
              if (i2 == i3) {
                continue;
              }
              if (i3 + 1 == exclude_atom) {
                continue;
              }
              if (0 == (ap3->color & HETERO_FLAG)) {
                continue;
              }
              for (int j = 2; j <= 6; j++) {
                if (0 == (length_matrix[i1][i2] & (1 << j))) {
                  continue;
                }
                for (int j1 = 2; j1 <= 5; j1++) {
                  if (0 == (length_matrix[i2][i3] & (1 << j1))) {
                    continue;
                  }
                  for (int j2 = 2; j2 <= 5; j2++) {
                    if (0 == (length_matrix[i3][i1] & (1 << j2))) {
                      continue;
                    }
                    // Make sure we have a real triangle
                    if (j + j1 == j2) {
                      continue;
                    }
                    if (j + j2 == j1) {
                      continue;
                    }
                    if (j1 + j1 == j) {
                      continue;
                    }
                    if (j + j1 + j2 > 11) {
                      continue;  // spider too large
                    }
                    if (j + j1 + j2 < 8) {
                      continue;  // spider too small
                    }
                    seed = CLASS_SPIDER_SEED;
                    seed = NEXT_SEED(seed, j + j1 + j2);
                    seed = NEXT_SEED(seed, j * j1 * j2);
                    // distinguish ring sizes of second ring-subst
                    // for (k=3; k<15; k++)
                    for (int k = 3; k < 9; k++) {
                      if (!(ap2->rsize_flags & (1 << k))) {
                        continue;
                      }
                      ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, k * 213));
                      result++;
                    }
                    // i1+1, ap1->rsize_flags, j,
                    // i2+1, ap2->rsize_flags, j1,
                    // i3+1, ap3->atom_symbol, j2);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  /* Collect bits that describe scaffolds. Those bits cannot (yet?) be used for
   * SSS screening */
  /* The method first collects a variant of the extended connectivity but
   * derived only        */
  /* on the ring bond graph. The bits are then set by combining atom features
   * with the        */
  /* corresponding scaffold extended connectivity of that atom. In addition, we
   * add           */
  /* bits to code connectivity of scaffolds by certain numbers of bonds using
   * the             */
  /* previously collected length_matrix. */
  if (!as_query && (which_bits & (USE_NON_SSS_BITS))) {
    seed = NON_SSS_SEED;
    extcon.assign(mp->n_atoms, 0);
    extcon2.assign(mp->n_atoms, 0);
    /* initialized extended connectivity */
    ap = mp->atom_array;
    for (int j = 0; j < mp->n_atoms; j++, ap++) {
      if (atom_status[j] <= 0) {
        continue;
      }
      extcon[j] = ap->rsize_flags;
    }
    /* propagate extended connectivity to neighbours for a few cycles */
    for (int i = 0; i < 32; i++) {
      for (int j = 0; j < mp->n_atoms; j++) {
        extcon2[j] = 0;
      }
      for (int j = 0; j < mp->n_atoms; j++) {
        /* skip non-ring atoms */
        if (atom_status[j] <= 0) {
          continue;
        }
        extcon2[j] = atom_status[j] * 3 + (extcon[j] * 0xF);
        sum = 0;
        prod = 0;
        for (int jj = 0; jj < nbp[j].n_ligands; jj++) {
          // only propagate through ring bonds
          if (bond_status[nbp[j].bonds[jj]] <= 0) {
            continue;
          }
          sum += extcon[nbp[j].atoms[jj]];
          prod *= 0xFF & (1 + extcon[nbp[j].atoms[jj]]);
        }
        extcon2[j] = ((extcon2[j] + (sum * 191) + prod) << 8) +
                     ((extcon2[j] & 0xFF0000) >> 16);
        extcon2[j] &= 0xFFFFFF;
      }
      for (int j = 0; j < mp->n_atoms; j++) {
        extcon[j] = extcon2[j];
      }
    }

    /* propagate smallest hash to all members of ring system */
    for (;;) {
      changed = FALSE;
      for (int j = 0; j < mp->n_atoms; j++) {
        /* skip non-ring atoms */
        if (atom_status[j] <= 0) {
          continue;
        }
        for (int jj = 0; jj < nbp[j].n_ligands; jj++) {
          // only propagate through ring bonds
          if (bond_status[nbp[j].bonds[jj]] <= 0) {
            continue;
          }
          if (extcon[j] < extcon[nbp[j].atoms[jj]]) {
            changed = TRUE;
            extcon[nbp[j].atoms[jj]] = extcon[j];
          }
        }
      }
      if (!changed) {
        break;
      }
    }

    /* Now, use extcon to set bits */
    ap = mp->atom_array;
    for (int j = 0; j < mp->n_atoms; j++, ap++) {
      if (j + 1 == exclude_atom) {
        continue;
      }
      /* skip non-ring atoms */
      if (atom_status[j] <= 0) {
        continue;
      }
      if (which_bits & USE_SCAFFOLD_IDS) {
        // ADD_BIT(fp_counts, ncounts, extcon[j]*1013 + seed);
        ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, extcon[j] * 10013));
        // ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, extcon[j]*40013));
        result += 1;
      }

      if (which_bits & USE_SCAFFOLD_LINKS) {
        /* needs to be a substitution point */
        if (degree[j] <= atom_status[j]) {
          continue;
        }
        for (int jj = 0; jj < mp->n_atoms; jj++) {
          if (jj + 1 == exclude_atom) {
            continue;
          }
          /* skip non-ring atoms */
          if (atom_status[jj] <= 0) {
            continue;
          }
          if (j == jj) {
            continue;
          }
          /* needs to be a substitution point */
          if (degree[jj] <= atom_status[jj]) {
            continue;
          }
          int pathLength = 12;
          for (int k = 0; k < 12; k++) {
            if (length_matrix[j][jj] == (1 << k)) {
              pathLength = k;
              break;
            }
          }
          if (pathLength >= 12) {
            continue;  // multiple paths => no chain connection
          }
          if (pathLength >= 3) {
            continue;  // only short one counts
          }
          // bit for ring system
          ADD_BIT(fp_counts, ncounts,
                  NEXT_SEED(NEXT_SEED(seed, 3 * pathLength),
                            extcon[j] * 1013 + extcon[jj] * 2003));
          // bit for position
          ADD_BIT(fp_counts, ncounts,
                  NEXT_SEED(NEXT_SEED(seed, 5 * pathLength),
                            extcon2[j] * 2013 + extcon[jj] * 1003));
          ADD_BIT(fp_counts, ncounts,
                  NEXT_SEED(NEXT_SEED(seed, 5 * pathLength),
                            extcon[j] * 2013 + extcon2[jj] * 1003));
          ADD_BIT(fp_counts, ncounts,
                  NEXT_SEED(NEXT_SEED(seed, 7 * pathLength),
                            extcon2[j] * 3013 + extcon2[jj] * 3003));
          result += 4;
        }
      }
    }
    if (which_bits & USE_SCAFFOLD_COLORS) {
      ap = mp->atom_array;
      for (int j = 0; j < mp->n_atoms; j++, ap++) {
        if (j + 1 == exclude_atom) {
          continue;
        }
        if (atom_status[j] <= 0) {
          continue;
        }
        tmp1 = 0;
        if (ap->atom_symbol == "C") {
          tmp1 = 101;
        } else if (ap->atom_symbol == "O") {
          tmp1 = 301;
        } else if (ap->atom_symbol == "N") {
          tmp1 = 401;
        } else if (ap->atom_symbol == "S") {
          tmp1 = 601;
        } else if (ap->atom_symbol == "P") {
          tmp1 = 701;
        } else if (AtomSymbolMatch(ap->atom_symbol, "F,Cl,Br,I,At")) {
          tmp1 = 901;
        }
        extcon[j] = ap->rsize_flags + tmp1;
      }
      for (int i = 0; i < 32; i++) {
        for (int j = 0; j < mp->n_atoms; j++) {
          extcon2[j] = 0;
        }
        for (int j = 0; j < mp->n_atoms; j++) {
          /* skip non-ring atoms */
          if (atom_status[j] <= 0) {
            continue;
          }
          extcon2[j] = atom_status[j] * 3 + (extcon[j] * 0xF);
          sum = 0;
          prod = 0;
          for (int jj = 0; jj < nbp[j].n_ligands; jj++) {
            // only propagate through ring bonds
            if (bond_status[nbp[j].bonds[jj]] <= 0) {
              continue;
            }
            sum += extcon[nbp[j].atoms[jj]];
            prod *= 0xFF & (1 + extcon[nbp[j].atoms[jj]]);
          }
          extcon2[j] = ((extcon2[j] + (sum * 191) + prod) << 8) +
                       ((extcon2[j] & 0xFF0000) >> 16);
          extcon2[j] &= 0xFFFFFF;
        }
        for (int j = 0; j < mp->n_atoms; j++) {
          extcon[j] = extcon2[j];
        }
      }
      /* propagate smallest hash to all members of ring system */
      for (;;) {
        changed = FALSE;
        for (int j = 0; j < mp->n_atoms; j++) {
          /* skip non-ring atoms */
          if (atom_status[j] <= 0) {
            continue;
          }
          for (int jj = 0; jj < nbp[j].n_ligands; jj++) {
            // only propagate through ring bonds
            if (bond_status[nbp[j].bonds[jj]] <= 0) {
              continue;
            }
            if (extcon[j] < extcon[nbp[j].atoms[jj]]) {
              changed = TRUE;
              extcon[nbp[j].atoms[jj]] = extcon[j];
            }
          }
        }
        if (!changed) {
          break;
        }
      }
      for (int j = 0; j < mp->n_atoms; j++) {
        if (j + 1 == exclude_atom) {
          continue;
        }
        ADD_BIT(fp_counts, ncounts,
                NEXT_SEED(NEXT_SEED(seed, 4), extcon[j] * 3013));
        // ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,4),
        // extcon[j]*7013)); ADD_BIT(fp_counts, ncounts, extcon[j]*16013 +
        // 4*seed);
        result += 1;
      }
    }
  }

  if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS)) {
    /* de-allocate length_matrix */
  }

  // seed = 0;
  // for (i=0; i<nbytes; i++)
  //    seed = (0xFF&fingerprint[i]) ^ ((seed>>8) | ((seed&0xFF)<<16));
  //         CountBits(fingerprint, nbytes), seed, result);
  // FORTIFY    fortifyTest(total_bytes_allocated,
  // "SetFingerprintCountsWithFocus");
  return (result);
}

#undef SINGLE
#undef DOUBLE
#undef TRIPLE
#undef AROMATIC
#undef NONE
#undef TRUE
#undef FALSE

// Enumerates the rings through ring atoms/bonds (up to max_size) and sets the
// ring size flags (bit k: member of a ring of size k; bit 0: any ring)
void markRingsRecursive(reaccs_molecule_t *mp, std::vector<int> &touchedAtoms,
                        std::vector<int> &touchedBonds, int startIndex,
                        int pathLength, int currentIndex, int maxSize,
                        const std::vector<neighbourhood_t> &nbp) {
  for (int i = 0; i < nbp[currentIndex].n_ligands; i++) {
    int ai = nbp[currentIndex].atoms[i];
    if (ai < startIndex) {
      continue;
    }
    if (ai == startIndex) {
      if (pathLength < 3) {
        continue;
      }
      for (int j = 0; j < mp->n_atoms; j++) {
        if (touchedAtoms[j]) {
          mp->atom_array[j].rsize_flags |= (1 << pathLength);
        }
      }
      for (int j = 0; j < mp->n_bonds; j++) {
        if (touchedBonds[j]) {
          mp->bond_array[j].rsize_flags |= (1 << pathLength);
        }
      }
      continue;
    }
    if (touchedAtoms[ai]) {
      continue;
    }
    if (pathLength + 1 > maxSize) {
      continue;
    }
    if (mp->atom_array[ai].rsize_flags == 0) {
      continue;
    }
    int bi = nbp[currentIndex].bonds[i];
    if (mp->bond_array[bi].rsize_flags == 0) {
      continue;
    }
    touchedAtoms[ai] = 1;
    touchedBonds[bi] = 1;
    markRingsRecursive(mp, touchedAtoms, touchedBonds, startIndex,
                       pathLength + 1, ai, maxSize, nbp);
    touchedAtoms[ai] = 0;
    touchedBonds[bi] = 0;
  }
}

}  // namespace

namespace {
using BondSet = std::vector<char>;

const int kSingle = 1;
const int kDouble = 2;
const int kTriple = 3;
const int kAromatic = 4;

// The ring sets used by the perception code: the (symmetrized) SSSR rings
// which only use usable bonds, followed by the rings obtainable by XORing two
// base rings which have exactly one connected path in common.
struct RingSets {
  std::vector<BondSet> base;
  std::vector<BondSet> pairs;
};

int cardinality(const BondSet &bs) {
  int res = 0;
  for (auto c : bs) {
    res += c != 0;
  }
  return res;
}

RingSets findRingSets(const ROMol &mol, const std::vector<char> &usableBond) {
  RingSets res;
  const auto nBonds = mol.getNumBonds();
  for (const auto &ring : mol.getRingInfo()->bondRings()) {
    bool ok = true;
    for (auto bi : ring) {
      ok = ok && usableBond[bi];
    }
    if (!ok) {
      continue;
    }
    BondSet bs(nBonds, 0);
    for (auto bi : ring) {
      bs[bi] = 1;
    }
    res.base.push_back(std::move(bs));
  }
  std::vector<int> touched(mol.getNumAtoms());
  for (size_t i = 0; i < res.base.size(); ++i) {
    for (size_t j = i + 1; j < res.base.size(); ++j) {
      std::fill(touched.begin(), touched.end(), 0);
      bool overlap = false;
      for (unsigned int b = 0; b < nBonds; ++b) {
        if (res.base[i][b] && res.base[j][b]) {
          overlap = true;
          const auto bond = mol.getBondWithIdx(b);
          ++touched[bond->getBeginAtomIdx()];
          ++touched[bond->getEndAtomIdx()];
        }
      }
      if (!overlap) {
        continue;
      }
      if (std::count(touched.begin(), touched.end(), 1) == 2) {
        BondSet bs(nBonds, 0);
        for (unsigned int b = 0; b < nBonds; ++b) {
          bs[b] = res.base[i][b] != res.base[j][b];
        }
        res.pairs.push_back(std::move(bs));
      }
    }
  }
  return res;
}

// port of PerceiveAromaticBonds(): six-rings (and fused rings with 4n+2
// bonds) which only consist of sp2 atoms are aromatic
void perceiveAromaticBonds(const ROMol &mol, std::vector<reaccs_bond_t> &bonds,
                           int nAtoms) {
  const int nBonds = bonds.size();
  std::vector<char> allBonds(nBonds, 1);
  auto rs = findRingSets(mol, allBonds);
  std::vector<char> inRing(nBonds, 0);
  for (const auto &r : rs.base) {
    for (int b = 0; b < nBonds; ++b) {
      inRing[b] |= r[b];
    }
  }
  std::vector<const BondSet *> rings;
  for (const auto &r : rs.pairs) {
    rings.push_back(&r);
  }
  for (const auto &r : rs.base) {
    rings.push_back(&r);
  }
  std::vector<int> spCount(nAtoms + 1);
  bool changed;
  do {
    changed = false;
    for (const auto ring : rings) {
      std::fill(spCount.begin(), spCount.end(), 0);
      int ndouble = 0, nsingle = 0;
      for (int b = 0; b < nBonds; ++b) {
        if ((*ring)[b] && bonds[b].bond_type == kSingle) {
          ++nsingle;
        }
      }
      bool isCumulene = false;
      for (int b = 0; b < nBonds; ++b) {
        if (!(*ring)[b]) {
          continue;
        }
        if (bonds[b].bond_type == kDouble) {
          ++ndouble;
          for (int k = 0; k < 2; ++k) {
            if (++spCount[bonds[b].atoms[k]] > 1) {
              isCumulene = true;
            }
          }
        } else if (bonds[b].bond_type == kTriple) {
          isCumulene = true;
        }
      }
      for (int b = 0; b < nBonds; ++b) {
        if ((*ring)[b] && bonds[b].bond_type == kAromatic) {
          for (int k = 0; k < 2; ++k) {
            if (spCount[bonds[b].atoms[k]] == 0) {
              ++spCount[bonds[b].atoms[k]];
            }
          }
        }
      }
      bool isAromatic = !isCumulene && ((cardinality(*ring) - 2) % 4) == 0;
      for (int b = 0; b < nBonds; ++b) {
        if ((*ring)[b] && (spCount[bonds[b].atoms[0]] != 1 ||
                           spCount[bonds[b].atoms[1]] != 1)) {
          isAromatic = false;
        }
      }
      if (isAromatic && (ndouble > 0 || nsingle > 0)) {
        for (int b = 0; b < nBonds; ++b) {
          if ((*ring)[b] && inRing[b] && bonds[b].bond_type != kAromatic) {
            changed = true;
            bonds[b].bond_type = kAromatic;
          }
        }
      }
    }
  } while (changed);
}

// port of PerceiveDYAromaticity(): Daylight-like aromaticity
void perceiveDYAromaticity(const ROMol &mol, std::vector<reaccs_atom_t> &atoms,
                           std::vector<reaccs_bond_t> &bonds,
                           const std::vector<neighbourhood_t> &nbp) {
  const int nAtoms = atoms.size();
  const int nBonds = bonds.size();
  std::vector<char> allBonds(nBonds, 1);
  auto allRings = findRingSets(mol, allBonds);
  if (allRings.base.empty()) {
    return;
  }
  std::vector<char> bondInRing(nBonds, 0);
  for (const auto &r : allRings.base) {
    for (int b = 0; b < nBonds; ++b) {
      bondInRing[b] |= r[b];
    }
  }
  std::vector<char> atomInRing(nAtoms, 0);
  std::vector<char> candidate(nAtoms, 0);
  for (int i = 0; i < nAtoms; ++i) {
    candidate[i] = atoms[i].atom_symbol != "C";
  }
  for (int b = 0; b < nBonds; ++b) {
    if (bonds[b].bond_type > kSingle && bonds[b].bond_type != kTriple) {
      candidate[bonds[b].atoms[0] - 1] = 1;
      candidate[bonds[b].atoms[1] - 1] = 1;
    }
  }
  std::vector<char> usable(nBonds, 0);
  for (int b = 0; b < nBonds; ++b) {
    if (!bondInRing[b]) {
      continue;
    }
    atomInRing[bonds[b].atoms[0] - 1] = 1;
    atomInRing[bonds[b].atoms[1] - 1] = 1;
    usable[b] =
        candidate[bonds[b].atoms[0] - 1] && candidate[bonds[b].atoms[1] - 1];
  }
  auto rs = findRingSets(mol, usable);
  if (rs.base.empty()) {
    return;
  }
  std::vector<BondSet> candidates;
  for (const auto &r : rs.base) {
    candidates.push_back(r);
  }
  for (const auto &r : rs.pairs) {
    candidates.push_back(r);
  }
  // fused pairs
  const size_t nRings = candidates.size();
  for (size_t i = 0; i < nRings; ++i) {
    for (size_t j = i + 1; j < nRings; ++j) {
      bool overlap = false;
      BondSet x(nBonds, 0);
      for (int b = 0; b < nBonds; ++b) {
        overlap = overlap || (candidates[i][b] && candidates[j][b]);
        x[b] = candidates[i][b] != candidates[j][b];
      }
      if (overlap && cardinality(x) == cardinality(candidates[i]) +
                                           cardinality(candidates[j]) - 2) {
        candidates.push_back(std::move(x));
      }
    }
  }

  const auto isSym = [&](int ai, const char *list) {
    return AtomSymbolMatch(atoms[ai].atom_symbol, list);
  };
  bool changed;
  do {
    changed = false;
    for (const auto &ring : candidates) {
      int npi = 0;
      bool conjugated = true;
      for (int i = 0; i < nAtoms; ++i) {
        if (!atomInRing[i]) {
          continue;
        }
        int inRingDouble = 0, inRingAromatic = 0;
        bool exoPull = false, isInRing = false;
        int localPi = 0;
        for (int j = 0; j < nbp[i].n_ligands; ++j) {
          const int bi = nbp[i].bonds[j];
          const int ai = nbp[i].atoms[j];
          if (bondInRing[bi] && ring[bi]) {
            isInRing = true;
            if (bonds[bi].bond_type == kAromatic) {
              ++inRingAromatic;
            }
            if (bonds[bi].bond_type == kDouble) {
              ++inRingDouble;
            }
          } else if (bonds[bi].bond_type == kDouble) {
            if (atoms[i].atom_symbol == "C" &&
                !bondInRing[bi] && isSym(ai, "O,S,P,N,L")) {
              exoPull = true;
            }
          }
        }
        if (!isInRing) {
          continue;
        }
        if ((inRingAromatic >= 1 || inRingDouble == 1) &&
            (isSym(i, "C,N,A,*") ||
             atoms[i].atom_symbol == "L")) {
          localPi = 1;
        } else if (inRingAromatic == 0 && inRingDouble == 0 &&
                   atoms[i].charge == 0 && isSym(i, "N,S,O")) {
          localPi = 2;
        } else if (inRingAromatic == 0 && inRingDouble == 0 &&
                   atoms[i].charge == 0 && exoPull &&
                   atoms[i].atom_symbol == "C") {
          localPi = 0;
        } else {
          conjugated = false;
        }
        if (atoms[i].charge < 0 && isSym(i, "C,N")) {
          conjugated = false;
        }
        npi += localPi;
      }
      if (!conjugated || npi % 4 != 2) {
        continue;
      }
      for (int b = 0; b < nBonds; ++b) {
        if (bondInRing[b] && ring[b] && bonds[b].bond_type != kAromatic) {
          bonds[b].bond_type = kAromatic;
          changed = true;
        }
      }
    }
  } while (changed);

  std::fill(candidate.begin(), candidate.end(), 0);
  for (const auto &b : bonds) {
    if (b.bond_type == kAromatic) {
      candidate[b.atoms[0] - 1] = candidate[b.atoms[1] - 1] = 1;
    }
  }
  for (int b = 0; b < nBonds; ++b) {
    if (candidate[bonds[b].atoms[0] - 1] && candidate[bonds[b].atoms[1] - 1] &&
        bondInRing[b] && bonds[b].bond_type == kSingle) {
      bonds[b].bond_type = kAromatic;
    }
  }
}

// port of the carbon part of GuessHCountsFromSubstitution() as it applies to
// molecules without substitution-count queries
void guessSubstitution(const ROMol &mol, std::vector<reaccs_atom_t> &atoms,
                       const std::vector<reaccs_bond_t> &bonds,
                       const std::vector<neighbourhood_t> &nbp) {
  for (size_t i = 0; i < atoms.size(); ++i) {
    if (atoms[i].charge != 0 ||
        mol.getAtomWithIdx(i)->getNumRadicalElectrons() != 0 ||
        atoms[i].atom_symbol != "C") {
      continue;
    }
    int nsingle = 0, ndouble = 0, ntriple = 0, naromatic = 0, nother = 0;
    for (int j = 0; j < nbp[i].n_ligands; ++j) {
      switch (bonds[nbp[i].bonds[j]].bond_type) {
        case kSingle:
          ++nsingle;
          break;
        case kDouble:
          ++ndouble;
          break;
        case kAromatic:
          ++naromatic;
          break;
        case kTriple:
          ++ntriple;
          break;
        default:
          ++nother;
      }
    }
    if (nother > 0) {
      continue;
    }
    if (naromatic == 2) {
      ++nsingle;
      ++ndouble;
      naromatic = 0;
    }
    if (naromatic > 0) {
      continue;
    }
    if (nsingle + 2 * ndouble + 3 * ntriple == 4) {
      atoms[i].sub_desc = SUB_AS_IS;
    }
  }
}
}  // namespace

std::vector<std::uint32_t> getAvalonCounts(const ROMol &mol,
                                           std::uint32_t fpSize,
                                           std::uint32_t bitFlags, bool isQuery,
                                           int focusAtom,
                                           bool accumulateAsQuery) {
  PRECONDITION(fpSize > 0, "fpSize must be > 0");
  PRECONDITION(
      focusAtom < 0 || static_cast<unsigned>(focusAtom) < mol.getNumAtoms(),
      "focusAtom out of range");
  std::vector<std::uint32_t> res(fpSize, 0);
  const int nAtoms = mol.getNumAtoms();
  if (!nAtoms) {
    return res;
  }

  // Avalon does its own aromaticity perception, so work with a kekulized copy
  auto tmol = std::make_unique<RWMol>(mol);
  try {
    MolOps::Kekulize(*tmol, true);
  } catch (const MolSanitizeException &) {
    // keep the aromatic bonds as they are
  }
  if (!tmol->getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(*tmol);
  }
  if (!isQuery) {
    for (const auto atom : tmol->atoms()) {
      if (atom->needsUpdatePropertyCache()) {
        tmol->updatePropertyCache(false);
        break;
      }
    }
  }
  const auto *lmol = tmol.get();
  const auto &ri = *lmol->getRingInfo();
  const int nBonds = lmol->getNumBonds();

  std::vector<reaccs_atom_t> atoms(nAtoms);
  std::vector<reaccs_bond_t> bonds(nBonds);
  std::vector<neighbourhood_t> nbp(nAtoms);
  std::vector<int> explicitH(nAtoms + 1, 0);
  std::vector<int> atomStatus(nAtoms, 0);
  std::vector<int> bondStatus(nBonds, 0);

  for (const auto atom : lmol->atoms()) {
    auto &a = atoms[atom->getIdx()];
    const auto anum = atom->getAtomicNum();
    if (anum == 0) {
      a.atom_symbol = "*";
    } else {
      a.atom_symbol =
          PeriodicTable::getTable()->getElementSymbol(anum);
    }
    a.charge = atom->getFormalCharge();
  }
  for (const auto bond : lmol->bonds()) {
    const int bi = bond->getIdx();
    const int bai = bond->getBeginAtomIdx();
    const int eai = bond->getEndAtomIdx();
    auto &b = bonds[bi];
    b.atoms[0] = bai + 1;
    b.atoms[1] = eai + 1;
    if (bond->getIsAromatic() || bond->getBondType() == Bond::AROMATIC) {
      b.bond_type = kAromatic;
    } else if (bond->getBondType() == Bond::SINGLE) {
      b.bond_type = kSingle;
    } else if (bond->getBondType() == Bond::DOUBLE) {
      b.bond_type = kDouble;
    } else if (bond->getBondType() == Bond::TRIPLE) {
      b.bond_type = kTriple;
    } else {
      b.bond_type = 8;
    }
    nbp[bai].atoms.push_back(eai);
    nbp[bai].bonds.push_back(bi);
    nbp[bai].n_ligands++;
    nbp[eai].atoms.push_back(bai);
    nbp[eai].bonds.push_back(bi);
    nbp[eai].n_ligands++;
    // explicit hydrogens in the graph count towards the H count
    if (lmol->getAtomWithIdx(eai)->getAtomicNum() == 1) {
      explicitH[bai + 1]++;
    } else if (lmol->getAtomWithIdx(bai)->getAtomicNum() == 1) {
      explicitH[eai + 1]++;
    }
    bondStatus[bi] = ri.numBondRings(bi);
    if (bondStatus[bi]) {
      atomStatus[bai]++;
      atomStatus[eai]++;
    }
  }
  const auto kekuleBonds = bonds;

  reaccs_molecule_t mp;
  mp.n_atoms = nAtoms;
  mp.n_bonds = nBonds;
  mp.atom_array = atoms.data();
  mp.bond_array = bonds.data();

  for (int i = 0; i < nAtoms; ++i) {
    atoms[i].rsize_flags = atomStatus[i] > 0 ? 1 : 0;
  }
  for (int i = 0; i < nBonds; ++i) {
    bonds[i].rsize_flags = bondStatus[i] > 0 ? 1 : 0;
  }
  std::vector<int> touchedAtoms(nAtoms, 0);
  std::vector<int> touchedBonds(nBonds, 0);
  for (int i = 0; i < nAtoms; ++i) {
    if (atoms[i].rsize_flags == 0) {
      continue;
    }
    touchedAtoms[i] = 1;
    markRingsRecursive(&mp, touchedAtoms, touchedBonds, i, 1, i, 14, nbp);
    touchedAtoms[i] = 0;
  }

  // The Avalon wrapper in AvalonTools.cpp runs the algorithm twice for
  // non-query molecules: once with Avalon's simple aromaticity model and once
  // with Daylight-like aromaticity, accumulating the results.
  std::vector<int> counts(fpSize, 0);
  const int nPasses = isQuery ? 1 : 2;
  for (int pass = 0; pass < nPasses; ++pass) {
    const bool dy = pass == 1;
    const bool queryMode = dy ? accumulateAsQuery : isQuery;
    for (int i = 0; i < nBonds; ++i) {
      bonds[i].bond_type = kekuleBonds[i].bond_type;
    }
    auto hCount = explicitH;
    for (auto &a : atoms) {
      a.sub_desc = 0;
    }
    if (!queryMode) {
      for (const auto atom : lmol->atoms()) {
        hCount[atom->getIdx() + 1] += atom->getTotalNumHs();
      }
    } else {
      guessSubstitution(*lmol, atoms, bonds, nbp);
    }
    if (dy) {
      perceiveDYAromaticity(*lmol, atoms, bonds, nbp);
    } else {
      perceiveAromaticBonds(*lmol, bonds, nAtoms);
    }
    CountFingerprintPatterns(&mp, nbp, hCount.data(), atomStatus.data(),
                             bondStatus.data(), counts.data(), fpSize, bitFlags,
                             queryMode, focusAtom + 1);
  }
  for (unsigned int i = 0; i < fpSize; ++i) {
    res[i] = counts[i] > 0 ? counts[i] : 0;
  }
  return res;
}

// ---- FingerprintGenerator integration ----

AvalonArguments::AvalonArguments(std::uint32_t bitFlags, bool isQuery,
                                 bool countSimulation,
                                 const std::vector<std::uint32_t> countBounds,
                                 std::uint32_t fpSize, bool accumulateAsQuery)
    : FingerprintArguments(countSimulation, countBounds, fpSize, 1, false),
      d_bitFlags(bitFlags),
      df_isQuery(isQuery),
      df_accumulateAsQuery(accumulateAsQuery) {
  PRECONDITION(fpSize > 0, "fpSize must be > 0");
}

std::string AvalonArguments::infoString() const {
  return "AvalonArguments bitFlags=" + std::to_string(d_bitFlags) +
         " isQuery=" + std::to_string(df_isQuery) +
         " accumulateAsQuery=" + std::to_string(df_accumulateAsQuery) + " " +
         commonArgumentsString();
}

void AvalonArguments::toJSON(boost::property_tree::ptree &pt) const {
  pt.put("type", "AvalonArguments");
  pt.put("bitFlags", d_bitFlags);
  pt.put("isQuery", df_isQuery);
  pt.put("accumulateAsQuery", df_accumulateAsQuery);
  FingerprintArguments::toJSON(pt);
}

void AvalonArguments::fromJSON(const boost::property_tree::ptree &pt) {
  d_bitFlags = pt.get<std::uint32_t>("bitFlags", d_bitFlags);
  df_isQuery = pt.get<bool>("isQuery", df_isQuery);
  df_accumulateAsQuery =
      pt.get<bool>("accumulateAsQuery", df_accumulateAsQuery);
  FingerprintArguments::fromJSON(pt);
}

template <typename OutputType>
OutputType AvalonAtomEnv<OutputType>::getBitId(
    FingerprintArguments *, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, AdditionalOutput *, const bool,
    const std::uint64_t) const {
  return d_bitId;
}

template <typename OutputType>
void AvalonAtomEnv<OutputType>::updateAdditionalOutput(AdditionalOutput *,
                                                       std::uint64_t) const {}

template <typename OutputType>
std::vector<AtomEnvironment<OutputType> *>
AvalonEnvGenerator<OutputType>::getEnvironments(
    const ROMol &mol, FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, int,
    const AdditionalOutput *, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, bool) const {
  PRECONDITION(arguments, "bad arguments");
  if ((fromAtoms && !fromAtoms->empty()) ||
      (ignoreAtoms && !ignoreAtoms->empty())) {
    throw ValueErrorException(
        "fromAtoms and ignoreAtoms are not supported by the Avalon "
        "fingerprint generator");
  }
  auto *args = dynamic_cast<AvalonArguments *>(arguments);
  PRECONDITION(args, "arguments must be AvalonArguments");

  // the bit ids are the slots of the (hashed) Avalon count vector; each slot
  // contributes one environment per count
  auto counts =
      getAvalonCounts(mol, args->d_fpSize, args->d_bitFlags, args->df_isQuery,
                      -1, args->df_accumulateAsQuery);
  std::vector<AtomEnvironment<OutputType> *> res;
  for (std::uint32_t i = 0; i < counts.size(); ++i) {
    for (std::uint32_t j = 0; j < counts[i]; ++j) {
      res.push_back(new AvalonAtomEnv<OutputType>(i));
    }
  }
  return res;
}

template <typename OutputType>
std::string AvalonEnvGenerator<OutputType>::infoString() const {
  return "AvalonEnvGenerator fpSize=" + std::to_string(d_fpSize);
}

template <typename OutputType>
void AvalonEnvGenerator<OutputType>::toJSON(
    boost::property_tree::ptree &pt) const {
  pt.put("type", "AvalonEnvGenerator");
  pt.put("fpSize", d_fpSize);
  AtomEnvironmentGenerator<OutputType>::toJSON(pt);
}

template <typename OutputType>
void AvalonEnvGenerator<OutputType>::fromJSON(
    const boost::property_tree::ptree &pt) {
  d_fpSize = pt.get<std::uint32_t>("fpSize", d_fpSize);
  AtomEnvironmentGenerator<OutputType>::fromJSON(pt);
}

template <typename OutputType>
OutputType AvalonEnvGenerator<OutputType>::getResultSize() const {
  return d_fpSize;
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getAvalonGenerator(
    const AvalonArguments &args) {
  return new FingerprintGenerator<OutputType>(
      new AvalonEnvGenerator<OutputType>(args.d_fpSize),
      new AvalonArguments(args), nullptr, nullptr, false, false);
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getAvalonGenerator(std::uint32_t fpSize,
                                                     std::uint32_t bitFlags,
                                                     bool isQuery,
                                                     bool accumulateAsQuery) {
  AvalonArguments args(bitFlags, isQuery, false, {1, 2, 4, 8}, fpSize,
                       accumulateAsQuery);
  return getAvalonGenerator<OutputType>(args);
}

template class AvalonAtomEnv<std::uint32_t>;
template class AvalonAtomEnv<std::uint64_t>;
template class AvalonEnvGenerator<std::uint32_t>;
template class AvalonEnvGenerator<std::uint64_t>;

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t> *
getAvalonGenerator(const AvalonArguments &);
template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *
getAvalonGenerator(const AvalonArguments &);
template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t> *
getAvalonGenerator(std::uint32_t, std::uint32_t, bool, bool);
template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *
getAvalonGenerator(std::uint32_t, std::uint32_t, bool, bool);

}  // namespace AvalonFP
}  // namespace RDKit
