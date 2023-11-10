#
# Copyright (C) 2007-2023 Greg Landrum and other RDKit contributors
#   All Rights Reserved
#
""" Implementation of the clustering algorithm published in:
  Butina JCICS 39 747-750 (1999)

"""
import numpy

from rdkit import RDLogger

logger = RDLogger.logger()


def EuclideanDist(pi, pj):
  dv = numpy.array(pi) - numpy.array(pj)
  return numpy.sqrt(dv * dv)


def _DenseClusterData(data, nPts, distThresh, isDistData, distFunc, reordering):
  if isDistData and len(data) > (nPts * (nPts - 1) / 2):
    logger.warning("Distance matrix is too long")

  # Create an empty matrix of booleans with a point for each i,j pair of molecule numbers.
  # Also tally the count of all pairs, skipping duplicates, and all self-comparison.
  # For efficiency, we never change the size of this bit bit-table.  We just flip bits.
  matrixPassingThreshold = numpy.zeros([nPts, nPts], dtype=bool)
  counts = numpy.zeros([nPts], dtype=numpy.int32)
  dmIdx = 0
  for i in range(nPts):
    for j in range(i):
      if not isDistData:
        dij = distFunc(data[i], data[j])
      else:
        dij = data[dmIdx]
        dmIdx += 1
      if dij <= distThresh:
        matrixPassingThreshold[i, j] = True
        matrixPassingThreshold[j, i] = True
        counts[i] += 1
        counts[j] += 1
      else:
        pass

  # sort by the number of neighbors:
  tLists = [(x, i) for i, x in enumerate(counts)]
  tLists.sort(reverse=True)

  res = []
  seen = numpy.zeros([nPts], dtype=bool)
  while tLists:
    _, idx = tLists.pop(0)
    if seen[idx]:
      continue

    newCluster = [idx]
    for jdxOtherCluster in range(nPts):
      if matrixPassingThreshold[idx, jdxOtherCluster]:
        if not seen[jdxOtherCluster]:
          newCluster.append(jdxOtherCluster)
          seen[jdxOtherCluster] = True

    # update the number of neighbors:
    # remove all members of the new cluster from the list of
    # neighbors and reorder the tLists
    if reordering:
      # get the list of affected molecules, i.e. all molecules
      # which have at least one of the members of the new cluster
      # as a neighbor
      neighbors = set()
      for idxNewCluster in newCluster:
        for jdxOtherCluster in range(nPts):
          if matrixPassingThreshold[idxNewCluster, jdxOtherCluster]:
            neighbors.add(jdxOtherCluster)
      neighbors = frozenset(neighbors)

      # loop over all remaining molecules in tLists but only
      # consider unassigned and affected compounds
      for tlistN, tlistCountAndIdx in enumerate(tLists):
        idxRemainingCluster = tlistCountAndIdx[1]
        if seen[idxRemainingCluster] or (idxRemainingCluster not in neighbors):
          continue
        # update the number of neighbors
        for jdxOtherCluster in newCluster:
          matrixPassingThreshold[idxRemainingCluster, jdxOtherCluster] = False
        tLists[tlistN] = (numpy.count_nonzero(matrixPassingThreshold[idxRemainingCluster, :]),
                          idxRemainingCluster)

      # now reorder the list
      tLists.sort(reverse=True)
    res.append(tuple(newCluster))
  return tuple(res)


def _SparseClusterData(data, nPts, distThresh, isDistData, distFunc, reordering):
  if isDistData and len(data) > (nPts * (nPts - 1) / 2):
    logger.warning("Distance matrix is too long")
  nbrLists = [None] * nPts
  for i in range(nPts):
    nbrLists[i] = []

  dmIdx = 0
  for i in range(nPts):
    for j in range(i):
      if not isDistData:
        dij = distFunc(data[i], data[j])
      else:
        dij = data[dmIdx]
        dmIdx += 1
      if dij <= distThresh:
        nbrLists[i].append(j)
        nbrLists[j].append(i)
  # sort by the number of neighbors:
  tLists = [(len(y), x) for x, y in enumerate(nbrLists)]
  tLists.sort(reverse=True)

  res = []
  seen = [0] * nPts
  while tLists:
    _, idx = tLists.pop(0)
    if seen[idx]:
      continue
    tRes = [idx]
    for nbr in nbrLists[idx]:
      if not seen[nbr]:
        tRes.append(nbr)
        seen[nbr] = 1
    # update the number of neighbors:
    # remove all members of the new cluster from the list of
    # neighbors and reorder the tLists
    if reordering:
      # get the list of affected molecules, i.e. all molecules
      # which have at least one of the members of the new cluster
      # as a neighbor
      nbrNbr = [nbrLists[t] for t in tRes]
      nbrNbr = frozenset().union(*nbrNbr)
      # loop over all remaining molecules in tLists but only
      # consider unassigned and affected compounds
      for x, y in enumerate(tLists):
        y1 = y[1]
        if seen[y1] or (y1 not in nbrNbr):
          continue
        # update the number of neighbors
        nbrLists[y1] = set(nbrLists[y1]).difference(tRes)
        tLists[x] = (len(nbrLists[y1]), y1)
      # now reorder the list
      tLists.sort(reverse=True)
    res.append(tuple(tRes))
  return tuple(res)


def ClusterData(data, nPts, distThresh, isDistData=False, distFunc=EuclideanDist, reordering=False,
                useSparseApproach=False):
  """  clusters the data points passed in and returns the list of clusters

    **Arguments**

      - data: a list of items with the input data (see discussion of
        _isDistData_ argument for the exception)

      - nPts: the number of points to be used

      - distThresh: elements within this range of each other are considered to
        be neighbors

      - isDistData: set this toggle when the data passed in is a
          distance matrix.  The distance matrix should be stored symmetrically.
          An example of how to do this:

            dists = [] for i in range(nPts):
              for j in range(i):
                dists.append( distfunc(i,j) )

      - distFunc: a function to calculate distances between points.
           Receives 2 points as arguments, should return a float

      - reordering: if this toggle is set, the number of neighbors is updated
           for the unassigned molecules after a new cluster is created such that
           always the molecule with the largest number of unassigned neighbors
           is selected as the next cluster center.
      
      - useSparseApproach: if this toggle is set, a sparse representation of the
        neighbor matrix will be used. This is slower, but uses much less memory,
        particularly for large numbers of points.


           
    **Returns**

      - a tuple of tuples containing information about the clusters:
         ( (cluster1_elem1, cluster1_elem2, ...),
           (cluster2_elem1, cluster2_elem2, ...), ...
         ) The first element for each cluster is its centroid.
  """
  if isDistData and len(data) > (nPts * (nPts - 1) / 2):
    logger.warning("Distance matrix is too long")

  if not useSparseApproach:
    return _DenseClusterData(data, nPts, distThresh, isDistData, distFunc, reordering)
  else:
    return _SparseClusterData(data, nPts, distThresh, isDistData, distFunc, reordering)
