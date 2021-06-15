"""
connect.py
Calculates the connections between poses of consecutive fragments based on overlap RMSD
During the calculation, poses of half of the fragments are clustered hierarchically, to speed up the calculation
Prints out the connectivity tree in JSON format

NOTE: First run get_msd_build.py to build the _get_msd Python extension module

Argument 1: nfrag, the number of fragments
Argument 2: the maximum RMSD
Argument 3: the maximum number of poses to consider (take the first poses in the .postatoms, .preatoms files)
Argument 4: the mimimum number of children clusters in 1st decomposition
Argument 5 to (4 + nfrag - 1): the "preatoms" portion of atom coordinates of the poses (the atoms that overlap with the next fragment)
Argument (5 + nfrag) to (5 + 2 * nfrag - 1): the "postatoms" portion of atom coordinates of the poses (the atoms that overlap with the previous fragment)
NOTE: the preatoms and postatoms are in .npy format, and must be sorted by ATTRACT rank! The first pose is rank 1, the 2nd is rank 2, etc.

Argument (4 + 2 * nfrag) - (4 + 3 * nfrag): optional: lists of pose indices to select for each fragment

Copyright 2015-2017 Sjoerd de Vries, Isaure Chauvot de Beauchene, TUM
"""

import sys, numpy as np, weakref
import json, itertools, bisect

#variables for computation parallelisation
MAXCHUNK = 4000000

#CLUSTERING = list of RMSD cutoffs for hierarchical clustering
#last CLUSTERING must be 0
#second-to-last is deredundant criterion (only one of two poses at this distance will be kept, as they are considered as identical)
CLUSTERING = [10, 8, 6, 5, 4, 3.5, 3, 2.5, 2, 1.7, 1.5, 1.0, 0.5, 0.1, 0]
MAX_CLUSTERING = len(CLUSTERING) - 1

# margin to consider if one member of a cluster could overlap with a member of another cluster
# based on the RMSD of the cluster representatives.
# 2 is the theoretical worse case; change to 9999 to disable all effects of clustering
CLUST_MARGIN = 2

###################################################
try:
  import cffi
  import _get_msd
  ffi = _get_msd.ffi
  def get_msd(c1, c2):
    #actually computes square-deviation (sd), not msd
    def npdata(a):
      return a.__array_interface__["data"][0]

    nc1 = len(c1)
    nc2 = len(c2)
    natom = c1.shape[1]
    assert c2.shape[1] == natom
    msd = np.empty((nc1, nc2), dtype=float)
    _get_msd.lib.get_msd(
      nc1, nc2, natom,
      ffi.cast("double *", npdata(c1)),
      ffi.cast("double *", npdata(c2)),
      ffi.cast("double *", npdata(msd)),
    )
    return msd
except ImportError:
  print >> sys.stderr, "Cannot find cffi, you will lose some speed"
  def get_msd(c1, c2):
    #actually computes square-deviation (sd), not msd
    d = c1[:, np.newaxis, :, :] - c2[np.newaxis, :, :, :]
    # dsq = d * d
    # msd = dsq.sum(axis=3).sum(axis=2)
    msd = np.einsum("...ijk,...ijk->...i", d,d)
    return msd

# Depth-first decomposition
def decompose(clusnr):
  best_conlevel = None # connection-level = priority of the cluster to be decomposed
                       # (here, correspond to depth)
  bestc = None # best cluster = cluster with highest priority
  clus = clusters[clusnr]
  for c in clus:
    if not len(c.children): continue
    if best_conlevel is None or c.conlevel > best_conlevel:
      best_conlevel = c.conlevel
      bestc = c
  if bestc is None: return False
  c = bestc
  clusters[clusnr].remove(c)
  if (clusnr % 2): # if preatoms
    c.decompose(fwd=True)        # check connections with postatoms of forward/downstream frag
    c.decompose_intra(fwd=False) # check connections with postatoms of same frag
  else: # if postatoms
    c.decompose_intra(fwd=True) # check connections with preatoms of same frag
    c.decompose(fwd=False)      # check connections with preatoms of backward/upstream frag
  for child in c.children:
    # additional ways to ensure depth search.
    conlevel = 0
    for con in itertools.chain(child.back_connections,child.connections):
      v = con.clusterlevel
      if v is not None and v > conlevel:
        conlevel = v
    child.conlevel = conlevel * child.clusterlevel
  clusters[clusnr].update(set(c.children))
  for cc in c.children:
    cc.check_deletion()
  return True

###################################################
class Cluster(object):
  __slots__ = ("clustid", "_splittable", "clusterlevel", "coors", "ranks", "all_ranks", "children", "nodes", "totstruc", "parent", "connections", \
    "back_connections", "_splitting", "_checking_delete", "conlevel", "_dead")
  def __init__(self, clustid, clusterlevel, coors, ranks):
    """
    A cluster can be in three forms:
    - leaf form. The cluster does not have any child clusters (yet).
        .coors contains all coordinates of the cluster (the first one is the representative).
        .ranks contains all ATTRACT ranks of the cluster.
    - node form. The cluster contains child clusters.
        .coors contains the coordinates of each representative of each child cluster.
        .ranks contains the ATTRACT rank of each representative of each child cluster.
    - singleton form. The cluster contains only a single member.
    There is a fourth form with coors = None, but this is not currently used in the code.
    The cluster is initially in leaf form (or in singleton form).
    """
    self._dead = False
    self.clustid = clustid # tuple
    # clustid of top-level clust of frag1-postatoms is (1,)
    #            2nd-level ........................ are (1,1), (1,2), ...
    #            top-level clust of frag1-preatoms  is (1001,)
    #            2nd-level ........................ are (1001,1), (1001,2)
    #            top-level ........ frag2-postatoms is (2,)
    #            top-level ........ frag2-preatoms  is (1002,)
    # ...
    self._splittable = True
    if coors is not None and len(coors) == 1: #singleton
      self.clusterlevel = MAX_CLUSTERING
    else:
      assert clusterlevel < MAX_CLUSTERING
      self.clusterlevel = clusterlevel #contains the clusterlevel of the cluster itself, not of the cluster children!
    if self.clusterlevel == MAX_CLUSTERING:
      self._splittable = False
    self.coors = coors #coordinates of the representative of the children clusters (or all coors if this clust is a leaf)
    self.ranks = ranks #ranks of the representative of the children clusters (or all ranks if this clust is a leaf)
    self.all_ranks = set(ranks) #ranks of the representative of the children clusters
    self.children = []
    self.nodes = 1 # nb of children
    if coors is not None:
      self.totstruc = coors.shape[0]
    self.parent = None
    self.connections = []
    self.back_connections = []
    self._splitting = False
    self._checking_delete = False
    self.conlevel = self.clusterlevel
    if self.clusterlevel is None:
      self.conlevel = -1
    if coors is not None and clusterlevel == MAX_CLUSTERING - 1: # deredundant level (cf line 17)
      r = self
      for cnr in range(len(self.coors)):
        c = Cluster(self.clustid + (cnr,), MAX_CLUSTERING, coors[cnr:cnr+1], ranks[cnr:cnr+1])
        c.parent = r
        self.children.append(c)
      self.nodes = len(self.children)
      self._splittable = False

  def _cluster(self, clusterlevel):
    '''subdivide cluster'''
    assert not self.children
    c = self.coors
    #clus: coordinates of the representative of each cluster
    clus = c[:1]
    #clus_indices: the coors indices of the structures of each cluster
    clus_indices = [[0]]
    chunksize = 20

    radius = CLUSTERING[clusterlevel]
    max_sd = radius * radius * c.shape[1]

    #This variable keeps track, for every structure in the chunk, into which new cluster it is sorted
    which_new_clust = np.zeros(chunksize, dtype=int)

    clustid = self.clustid
    if clustid is None: clustid = ()
    for n in range(1, len(c), chunksize):
      chunk = c[n:n+chunksize]

      #intra-chunk msd
      d = chunk[:, np.newaxis, :, :] - chunk[np.newaxis, :, :, :]
      intra_msd = np.einsum("...ijk,...ijk->...i", d,d)

      #chunk-cluster msd: compare to existing reference structures
      d = chunk[:, np.newaxis, :, :] - clus[np.newaxis, :, :, :]
      inter_msd = np.einsum("...ijk,...ijk->...i", d,d)

      for nn in range(len(chunk)):
        sort_clust = None
        # identify sets of structures that should go into the same new cluster
        # if they do not go in an existing cluster
        # => avoid creating to many references
        close_intra_clusts = (intra_msd[nn] < max_sd).nonzero()[0]
        intra_new_clusts = [which_new_clust[k] for k in close_intra_clusts if k < nn and which_new_clust[k] != -1]
        if len(intra_new_clusts):
          sort_clust = min(intra_new_clusts)
        close_inter_clusts = (inter_msd[nn] < max_sd).nonzero()[0]
        if len(close_inter_clusts):
          sort_clust2 = min(close_inter_clusts)
          if sort_clust is None or sort_clust > sort_clust2:
            sort_clust = sort_clust2
        if sort_clust is None:
          #new cluster
          which_new_clust[nn] = len(clus)
          clus = np.append(clus, chunk[nn][np.newaxis,:,:], axis=0)
          clus_indices.append([n+nn])
        else:
          clus_indices[sort_clust].append(n+nn)
          which_new_clust[nn] = -1

    indices = [i[0] for i in clus_indices]

    # Re-cluster to the closest reference
    clus_indices = [i[:1] for i in clus_indices]
    for n in range(0, len(c), chunksize):
      chunk = c[n:n+chunksize]
      d = chunk[:, np.newaxis, :, :] - clus[np.newaxis, :, :, :]
      inter_msd = np.einsum("...ijk,...ijk->...i", d,d)
      sort_clusts = np.argmin(inter_msd, axis=1)
      for nn in range(len(chunk)):
        if (n+nn) in indices: continue
        sort_clust = sort_clusts[nn]
        clus_indices[sort_clust].append(n+nn)

    for cnr,c in enumerate(clus_indices):
      ind = clus_indices[cnr]
      coors = self.coors[ind]
      ranks = self.ranks[ind]
      c = Cluster(clustid+(cnr+1,), clusterlevel, coors, ranks)
      self.children.append(c)
    self.nodes = len(self.children)
    self.totstruc = sum([c.totstruc for c in self.children])
    return clus, indices

  def cluster(self, clusterlevel):
    """Converts a cluster in leaf form to a cluster in node form."""
    assert clusterlevel < MAX_CLUSTERING
    clus, indices = self._cluster(clusterlevel)
    self.coors = clus
    self.ranks = self.ranks[indices]
    r = self
    for c in self.children:
      c.parent = r
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False

  def dissolve(self, clusterlevel):
    """Dissolve all direct child clusters, clustering them and linking all of their children to us"""
    self.ranks = np.concatenate([c.ranks for c in self.children])
    newchildren = []
    coors = []
    while len(self.children):
      child = self.children.pop()
      child.cluster(clusterlevel)
      coors.append(child.coors)
      newchildren += child.children
      self.totstruc = child.totstruc
    self.coors = np.concatenate(coors, axis=0)
    self.children = newchildren
    r = self
    for c in self.children:
      c.parent = r
    self.nodes = len(newchildren)
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False

  def reorganize(self):
    """Prune superfluous levels of clustering by dissolving child clusters with less than MINCHILD children
    Evokes reorganize() also on our children """
    for c in self.children:
      c.reorganize()
    if not len(self.children): return
    if len(self.children) >= MINCHILD: return
    oldchildren = [c for c in self.children if len(c.children)]
    if not len(oldchildren): return
    while len(self.children) < MINCHILD and len(oldchildren):
      child = oldchildren.pop(0)
      self.children.remove(child)
      self.children += child.children
      oldchildren += [c for c in child.children if len(c.children)]
    r = self
    coors = [c.coors[0] for c in self.children]
    self.coors = np.array(coors)
    for c in self.children:
      c.parent = r
    self.nodes = sum([c.nodes for c in self.children])
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False

  def split(self):
    """Sub-cluster ourselves and then our children, until we are a singleton or have more than 1 child
    Returns whether or not we are still splittable further"""
    if not self.children:
      assert self.clusterlevel is not None
      if self.clusterlevel == MAX_CLUSTERING: return False
      self.clusterlevel += 1
      if self.clusterlevel == MAX_CLUSTERING:
        self._splittable = False
        return False
      self.cluster(self.clusterlevel)
      r = self
      while len(self.children) == 1:
        child = self.children[0]
        self.clusterlevel = child.clusterlevel
        if self.clusterlevel >= MAX_CLUSTERING - 1:
          self._splittable = False
          break
        self.clusterlevel += 1
        child.cluster(self.clusterlevel)
        children = child.children
        for c in children: c.parent = r
        self.coors = child.coors
        self.ranks = child.ranks
        self.children = children
        self.nodes = len(children)
        self.totstruc = sum([c.totstruc for c in children])

      self.parent.add_nodes(self.nodes - 1)
      for c in self.children:
        if c._splittable:
          break
      else:
        self._splittable = False
      return True
    else:
      self._splitting = True
      oldnodes = self.nodes
      ok = False
      for c in self.children:
        c_splittable = c._splittable
        has_split = c.split()
        if has_split: ok = True
      newnodes = self.nodes
      self._splitting = False
      if self.parent is not None and newnodes > oldnodes:
        self.parent.add_nodes(newnodes-oldnodes)
      for c in self.children:
        if c._splittable:
          break
      else:
        self._splittable = False
      return ok

  def add_nodes(self, nodes):
    self.nodes += nodes
    if self.parent is not None and not self._splitting:
      self.parent.add_nodes(nodes)

  def check_deletion(self):
    """Check if cluster is dead because of lack of connections
    If so, remove the cluster and propagate check_deletion along the tree
    """
    if self._dead: return
    if self._checking_delete: return
    self._checking_delete = True
    while 1:
      has_c1, has_c2 = len(self.back_connections), len(self.connections)
      if has_c1 and has_c2: break # Forward and backward connections, cluster is alive

      # If not, check position in the chain, to know if we need
      # both forward and backward connections
      frag0 = self.clustid[0]
      if frag0 > 1000:
        pos = 2 * (frag0-1001) + 1
      else:
        pos = 2 * (frag0-1)
      if pos == 0 and has_c2: break # we are at the beginning of chain
                                    # and we do have fwd connections
      if pos == len(clusters) - 1 and has_c1: break # we are at the end of chain
                                                    # and we do have bwd connections

      #Cluster is dead. Remove it, remove all connections, and invoke check_deletion on the connected clusters
      self._dead = True
      clusters[pos].remove(self)
      # remove all our connection,
      # check if the thereby disconnected clusters are dead
      for o in list(self.back_connections):
        o.connections.remove(self)
        o.check_deletion()
      for o in list(self.connections):
        o.back_connections.remove(self)
        o.check_deletion()
      break
    self._checking_delete = False

  def decompose(self, fwd):
    '''
    When a node-form parent cluster is subdivided into children clusters,
    propagate from parent to children the connections with downstream(fwd=True) or upstream(fwd=False) clusters that are in overlap range
    then delete parent and connect children to the grand-parent (which is normally root)
    '''
    c1 = self.coors   # coordinates of the representatives of the children clusters at next level.
    if fwd:
      if not self.connections: return
      others = list(self.connections)
      # c2 = coord of the representative of the connected clusters of fwd frag
      c2 = np.concatenate([c.coors[0] for c in self.connections]).reshape(len(self.connections), len(self.connections[0].coors[0]), 3)
    else:
      if not self.back_connections: return
      others = list(self.back_connections)
      # c2 = coord of the representative of the connected clusters of bwd frag
      c2 = np.concatenate([c.coors[0] for c in self.back_connections]).reshape(len(self.back_connections), len(self.back_connections[0].coors[0]), 3)

    # For each of the representative of the connected clusters of fwd/bwd fragment (=c2),
    # check connectivity with each of the representative of the children clusters (= c1)
    # Divide c2 into chunks for memory saving
    chunksize = MAXCHUNK/len(c1)
    for chunkpos in range(0, len(c2), chunksize):
      # Account for clustering cutoff to check if any pose in cluster1 could
      # overlap with any pose in cluster2
      c_max_rmsd = [CLUSTERING[child.clusterlevel] * CLUST_MARGIN for child in self.children]
      max_rmsd0 = np.fromiter(c_max_rmsd,count=len(c_max_rmsd),dtype=float)[:, np.newaxis]
      max_rmsd2 = max_rmsd0 + max_rmsd # max_rmsd = cutoff given by user
      max_sd = (max_rmsd2**2) * c1.shape[1]

      c2_chunk = c2[chunkpos:chunkpos+chunksize]
      others_chunk = others[chunkpos:chunkpos+chunksize]
      if fwd:
        ocon = [o.back_connections for o in others_chunk]
        childcon = [c.connections for c in self.children]
      else:
        ocon = [o.connections for o in others_chunk]
        childcon = [c.back_connections for c in self.children]

      for o in ocon:
        o.remove(self)

      msd = get_msd(c1, c2_chunk)
      ch = self.children
      old_childnr = None

      for childnr, onr in zip(*np.where(msd < max_sd)):
        if childnr != old_childnr:
          c_child = ch[childnr]
          c_childcon = childcon[childnr]
          old_childnr = childnr
        ocon[onr].append(c_child)
        c_childcon.append(others_chunk[onr])

      for o in others:
        o.check_deletion()

  def decompose_intra(self, fwd):
    if fwd:
      if not self.connections: return
      others = list(self.connections)
    else:
      if not self.back_connections: return
      others = list(self.back_connections)

    if fwd:
      ocon = [o.back_connections for o in others]
      childcon = [c.connections for c in self.children]
    else:
      ocon = [o.connections for o in others]
      childcon = [c.back_connections for c in self.children]
    for o in ocon:
      o.remove(self)

    for childnr, child in enumerate(self.children):
      cranks = child.all_ranks
      for onr, o in enumerate(others):
        if cranks.intersection(o.all_ranks):
          ocon[onr].append(child)
          childcon[childnr].append(o)

    for o in others:
      o.check_deletion()

  def verify(self):
    ''' check that all kept connections have indeed low overlap rmsd'''
    if len(self.children): return
    if self.totstruc > 1: return
    cons = [con for con in self.connections if not len(con.children) and con.totstruc == 1]
    if not len(cons): return
    concoors = np.concatenate([con.coors[:1] for con in cons])
    c2 = self.coors[0]
    c1 = concoors
    d = c1 - c2
    msd = np.einsum("...ijk,ijk->...i", d,d)
    msd_low = (msd<(max_rmsd**2*c2.shape[0]))
    for n in range(len(cons)):
      assert msd_low[n]
    return

  def all_children(self):
    if not len(self.children):
      yield self
    else:
      for cc in self.children:
        for v in cc.all_children():
          yield v
###################################################
if __name__ == "__main__":
  nfrags = int(sys.argv[1])
  max_rmsd = float(sys.argv[2]) #overlapping cutoff (<3.0A recommended)
  maxstruc = int(sys.argv[3]) #take only the maxstruc top-ranked poses.
  MINCHILD = int(sys.argv[4])
  nargs = len(sys.argv) - 5
  # you can give selections of poses as arguments. see l38
  assert nargs == 2 * nfrags or nargs == 3 * nfrags, (nargs, nfrags)
  print >> sys.stderr, "NFRAGS", nfrags
  assert nfrags >= 2
  preatoms = sys.argv[5:nfrags+5]
  postatoms = sys.argv[nfrags+5:2*nfrags+5]
  print >> sys.stderr, "PREATOMS", preatoms
  print >> sys.stderr, "POSTATOMS", postatoms
  postatoms = [np.load(f) for f in postatoms]
  preatoms = [np.load(f) for f in preatoms]

  # lists of ranks to consider for each pose pool, counting from 1
  selections = [[] for n in range(nfrags)]
  if nargs == 3 * nfrags:
    selections = sys.argv[2*nfrags+5:3*nfrags+5]
    print >> sys.stderr, "SELECTIONS", selections
    selections = [np.array(sorted([int(l) for l in open(f) if len(l.strip())])) for f in selections]

  # dimensions = (Nb poses, Nb atoms * 3 coordinates)
  '''
  for a in postatoms: 
     assert len(a.shape) == 2 and a.shape[1] % 3 == 0, a.shape
  for a in preatoms:
    assert len(a.shape) == 2 and a.shape[1] % 3 == 0, a.shape

  for arr in (preatoms, postatoms):
    for anr, a in enumerate(arr):
      ncoor = a.shape[1] / 3
      arr[anr] = a.reshape(a.shape[0], ncoor, 3)
  '''

  #If you use both maxstruc and selection,
  #remove from selection what is beyond rank maxstruc
  if maxstruc > 0:
    for arr in (preatoms, postatoms):
      for anr, a in enumerate(arr):
        arr[anr] = arr[anr][:maxstruc]
    for selnr, sel in enumerate(selections):
      if not len(sel): continue
      pos = bisect.bisect_right(sel, maxstruc)
      selections[selnr] = sel[:pos]
      assert len(selections[selnr])

  # preatoms_frag(i) and postatoms_frag(i) must have
  # the same number of poses (nstruc)
  nstruc = []
  for n in range(nfrags):
    a1 = preatoms[n]
    a2 = postatoms[n]
    assert a1.shape[0] == a2.shape[0], (n, a1.shape, a2.shape)
    nstruc.append(a1.shape[0])

  # postatoms_frag(i) and preatoms_frag(i-1) must have
  # the same number of atoms, as they overlap in sequence.
  for n in range(1, nfrags):
    a1 = preatoms[n-1]
    a2 = postatoms[n]
    assert a1.shape[1] == a2.shape[1], (n, a1.shape, n+1, a2.shape)

  # Check that the pose exists for each rank in selection.
  for selnr, sel in enumerate(selections):
    for conf in sel:
      assert conf > 0 and conf <= nstruc[selnr], (conf, nstruc[selnr])

  ranks = [np.arange(s)+1 for s in nstruc]
  for n in range(nfrags):
    if len(selections[n]):
      preatoms[n] = preatoms[n][selections[n]-1]
      postatoms[n] = postatoms[n][selections[n]-1]
      ranks[n] = selections[n]
      nstruc[n] = len(selections[n])

  #Build cluster tree
  clusters = []
  for n in range(nfrags):
    for atoms in (postatoms, preatoms):
      i = n + 1
      if atoms is preatoms: i += 1000 # see clustid scheme (line 112)
      c = Cluster((i,), None, atoms[n], ranks[n])
      if not (n % 2):
        a = atoms[n]
        r = ranks[n]
        for nn in range(len(a)):
          cc = Cluster((i,nn), MAX_CLUSTERING, a[nn:nn+1], np.array(r[nn:nn+1]))
          c.children.append(cc)
        c.nodes = len(a)
        clusters.append([c])
        continue

      clusterlevel = 0
      c.cluster(clusterlevel)
      for clusterlevel in range(1, len(CLUSTERING)-1):
        if len(c.children) >= MINCHILD: break
        c.dissolve(clusterlevel)
      count = 0
      assert c.clusterlevel is None

      def split_all(c):
        """Split c and all of its children, all the way down"""
        global count
        if not c._splittable: return
        if not len(c.children):
          ok = c.split()
          count += 1
          if not ok: return
        for cc in c.children:
          split_all(cc)

      split_all(c)
      c.reorganize()
      print >> sys.stderr, n+1, nstruc[n], CLUSTERING[clusterlevel], len(c.children), c.nodes
      assert c.clusterlevel is None
      clusters.append([c])

  #Initialize tree connections, intra-fragment, flat
  for n in range(0, 2 * nfrags, 4):
    c1, c2 = clusters[n][0], clusters[n+1][0]
    #print >> sys.stderr,  n, n+1, len(c1.children), len(c2.children)
    for nn in range(len(c1.children)):
      cc1, cc2 = c1.children[nn], c2.children[nn]
      cc1.connections = [cc2]
      cc2.back_connections = [cc1]

  #Initialize tree connections, intra-fragment, hierarchical
  for n in range(2, 2 * nfrags, 4):
    c1, c2 = clusters[n][0], clusters[n+1][0]
    c1.connections.append(c2)
    c2.back_connections.append(c1)

  #Initialize tree connections, inter-fragment, flat to hierarchical
  for n in range(2, 2 * nfrags, 4):
    c1, c2 = clusters[n-1][0], clusters[n][0]
    for cc in c1.children:
      cc.connections = [c2]
      c2.back_connections.append(cc)
    clusters[n-1] = c1.children

  #Initialize tree connections, inter-fragment, hierarchical to flat
  for n in range(4, 2 * nfrags, 4):
    c1, c2 = clusters[n-1][0], clusters[n][0]
    for cc in c2.children:
      cc.back_connections = [c1]
      c1.connections.append(cc)
    clusters[n] = c2.children

  clusters[0] = clusters[0][0].children
  if len(clusters[-2]) > 1:
    clusters[-1] = clusters[-1][0].children

  clusters = [set(c) for c in clusters]


  #Decompose tree (divide clusters in sub-clusters)
  # Decompose first the clusters of poses for the extrem fragments
  # (fisrt and last in chain, then 2nd first and 2nd last...)
  # to gain efficiency, as less poses will be connected for those fragments
  # (at least in the exemples tested so far)
  step = 0
  to_decompose = []
  to_decompose0 = list(range(2, len(clusters), 4))
  while len(to_decompose0):
    to_decompose.append(to_decompose0.pop(0))
    if not len(to_decompose0): break
    to_decompose.append(to_decompose0.pop(-1))

  # Decompose preatoms and postatoms alternately at each clustering level,
  # so that you eliminate clusters that have no connections
  # e.g. from the postatoms to the preatoms, because the corresponding preatoms
  # had been eliminated as they did not connect to any pose of the downstream fragment
  for clusnr in to_decompose:
    done1, done2 = False, False
    while 1:
      if not (step % 5):
        print >> sys.stderr, [len(c) for c in clusters]
        print >> sys.stderr, [sum([cc.nodes for cc in c]) for c in clusters]
      step += 1
      if not done1:
        ok1 = decompose(clusnr)
        if not ok1:
          if done2: break
          done1 = True
      if not done2:
        ok2 = decompose(clusnr+1)
        if not ok2:
          if done1: break
          done2 = True

  print >> sys.stderr, [len(c) for c in clusters[::2]]
  print >> sys.stderr, [sum([cc.nodes for cc in c]) for c in clusters[::2]]

  #Verification
  for c in clusters[1::2]:
    for cc in c:
      cc.verify()

  #Sort clusters
  clusters = [list(c) for c in clusters]
  for c in clusters:
    c.sort(key=lambda clus: clus.ranks[0])

  #write out tree
  tree = {"nfrags":nfrags, "max_rmsd": max_rmsd, "clusters": [], "interactions": []}
  for cnr in range(1,len(clusters), 2):
    c = clusters[cnr]
    clus = []
    inter = []
    if cnr < len(clusters)-1:
      nex = clusters[cnr+1]
      nexmap = {}
      for vnr, v in enumerate(nex):
        nexmap[v] = vnr
    for ccnr, cc in enumerate(c):
      cclus = {"radius":CLUSTERING[cc.clusterlevel], "ranks": cc.ranks.tolist()}
      clus.append(cclus)
      if cnr < len(clusters)-1:
        for other in cc.connections:
          index = nexmap[other]
          inter.append((ccnr, index))
    inter.sort(key=lambda i: 100000*i[0]+i[1])
    tree["clusters"].append(clus)
    if cnr < len(clusters)-1:
      tree["interactions"].append(inter)

  json.dump(tree, sys.stdout)
