"""=================================================================================================
Read the distances calculated by fingerprint_distance. Expected format
fpt1    fpt2    Jaccard Bray-Curtis

First cluster in to connected components, then run UPGMA on each cluster

Michael Gribskov     06 May 2022
================================================================================================="""
import sys


def read_distance(filename):
    """---------------------------------------------------------------------------------------------
    Read distances from fingerprint_distance.py.  Distances are stored in a list of dict
    distance = [ {fpt1, fpt2, jaccard, bray_curtis}, ...]
    
    :param filename: string, name of file with distances
    :return: distance (list of dict), maximum (float), minimum (float)
    ---------------------------------------------------------------------------------------------"""
    try:
        distance = open(filename, 'r')
    except OSError:
        sys.stderr.write(f'distance_cluster - read_distance: cannot open distance file{filename}')

    maximum = {'jaccard': 0, 'bray-curtis': 0}
    minimum = {'jaccard': 1000000, 'bray-curtis': 1000000}

    distance_n = 0
    distance_list = []
    for line in distance:
        if line.startswith('#'):
            continue

        distance_n += 1
        fpt1, fpt2, jaccard, bray_curtis = line.rstrip().split('\t')
        maximum['jaccard'] = max(maximum['jaccard'], float(jaccard))
        minimum['jaccard'] = min(minimum['jaccard'], float(jaccard))
        maximum['bray-curtis'] = max(maximum['bray-curtis'], float(bray_curtis))
        minimum['bray-curtis'] = min(minimum['bray-curtis'], float(bray_curtis))
        distance_list.append({'fpt1':        fpt1,
                              'fpt2':        fpt2,
                              'jaccard':     float(jaccard),
                              'bray-curtis': float(bray_curtis)})

    return distance_list, maximum, minimum


def connected(distance, threshold, dtype='jaccard'):
    """---------------------------------------------------------------------------------------------
    Read distances from fingerprint_distance.py

    :param distance: list of dict, keys: fpt1, fpt2, jaccard, bray-curtis
    :return: list of string, IDs of fingerprints in connected components
    ---------------------------------------------------------------------------------------------"""
    cluster = []
    index = {}
    for d in distance:
        c1 = cluster_get(cluster, index, d['fpt1'])
        c2 = cluster_get(cluster, index, d['fpt2'])
        if d[dtype] > threshold:
            cluster_merge(cluster, index, c1, c2)

    return cluster, index


def cluster_get(cluster, index, id):
    """---------------------------------------------------------------------------------------------
    Return the cluster to which id belongs, if id is new, create a new cluster

    :param cluster: list of ids, clusters
    :param index: dictionary, key = id, cluster to which id belongs
    :param id: string, id of a fpt
    :return: int, cluster index
    ---------------------------------------------------------------------------------------------"""
    cid = None

    if id in index:
        cid = index[id]
    else:
        cid = len(cluster)
        index[id] = cid
        cluster.append([id])

    return cid


def cluster_merge(cluster, index, c1, c2):
    """---------------------------------------------------------------------------------------------
    merge cluster 1 and cluster 2, c1 and c2 are indices in cluster

    :param cluster: list of ids, clusters
    :param index: dictionary, key = id, cluster to which id belongs
    :param c1: int, index of cluster 1 in cluster
    :param c2: int, index of cluster 2 in cluster
    :return:
    ---------------------------------------------------------------------------------------------"""
    if c1 == c2:
        return

    for member in cluster[c2]:
        index[member] = c1

    cluster[c1] += cluster[c2]
    cluster[c2] = []

    return


def distance_matrix(distance, cluster, dtype='jaccard'):
    """---------------------------------------------------------------------------------------------
    construct a distance matrix for this cluster. Distances are store in a list of dict
    distance = [ {fpt1, fpt2, jaccard, bray_curtis}, ...]. The distance matrix is ordered the same
    as the cluster, and distances are assumed to be symmetric.

    :param distance: list of dict, keys: fpt1, fpt2, jaccard, bray-curtis
    :param cluster: list of string, IDs of fingerprints
    :return: list of list of float, distance matrix
    ---------------------------------------------------------------------------------------------"""
    dmat = [[0 for i in range(len(cluster))] for j in range(len(cluster))]
    for d in range(len(distance)):
        f1 = distance[d]['fpt1']
        f2 = distance[d]['fpt2']
        if f1 in cluster and f2 in cluster:
            i = cluster.index(f1)
            j = cluster.index(f2)
            dmat[i][j] = dmat[j][i] = distance[d][dtype]

    return dmat


class Tree:
    """=============================================================================================

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.idx = None
        self.id = ''
        children = []


class Upgma:
    """=============================================================================================

    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.dtype = 'jaccard'
        self.dmat = []
        self.tree = []

    def active(self):
        """-----------------------------------------------------------------------------------------
        Return a list of the active indices in self.tree (which are also the active rows and columns
        in self.dmat

        :return:
        -----------------------------------------------------------------------------------------"""
        idx = []
        for node in self.tree:
            idx.append(node.idx)

        return idx

    def load(self, distance, cluster):
        """-----------------------------------------------------------------------------------------
        Load the distance matrix and initialize tree leaves

        :param distance: list of dict, keys: fpt1, fpt2, jaccard, bray-curtis
        :param cluster: list of string, IDs fo fingerprints
        :return: int, number of taxa
        -----------------------------------------------------------------------------------------"""
        tree = self.tree
        taxa_n = 0
        for id in cluster:
            leaf = Tree()
            leaf.id = id
            leaf.idx = taxa_n
            tree.append(leaf)
            taxa_n += 1

        self.dmat = [[0 for i in range(len(cluster))] for j in range(len(cluster))]
        dmat = self.dmat
        dtype = self.dtype

        for d in range(len(distance)):
            f1 = distance[d]['fpt1']
            f2 = distance[d]['fpt2']
            if f1 in cluster and f2 in cluster:
                i = cluster.index(f1)
                j = cluster.index(f2)
                dmat[i][j] = dmat[j][i] = distance[d][dtype]

        return taxa_n

    def dmat_format(self, decimal=2, space=2):
        """-----------------------------------------------------------------------------------------
        Return a string with the formatted distance matrix

        :param decimal: int, number of digits past the decimal point
        :return: string, formatted matrix
        -----------------------------------------------------------------------------------------"""
        dmat = self.dmat
        dmatstr = ''
        formatstr = f'.{decimal}f'
        size = len(dmat)
        width = 0
        for i in range(size):
            for j in range(i, size):
                width = max(width, len(f'{dmat[i][j]:{formatstr}}'))
        # print(f'field width={width}')

        formatstr = f'{width + space}.{decimal}f'
        for i in range(size):
            dmatstr += ' ' * i * (width + space)
            for j in range(i, size):
                dmatstr += f'{dmat[i][j]:{formatstr}}'
            dmatstr += '\n'

        return dmatstr

    def smallest(self, initval=10):
        """---------------------------------------------------------------------------------------------
        find the smallest value in the current distance matrix, only look at values in groups, other
        indices have been merged

        :return: int, int, indices of the smallest distance
        ---------------------------------------------------------------------------------------------"""
        smallest = initval
        dmat = self.dmat
        tree = self.tree
        active = self.active()

        srow = 0
        scol = 0
        for i in range(len(active)):
            row = active[i]
            for j in range(i + 1, len(active)):
                col = active[j]
                if dmat[row][col] < smallest:
                    smallest = dmat[row][col]
                    srow = row
                    scol = col

        print(smallest, srow, scol)
        return srow, scol

    def mergetaxa(self, row, col):
        """-----------------------------------------------------------------------------------------
        merge the row and column taxa with the smallest value

        :param row: index of row (will be kept)
        :param col: index of column (will be merged)
        :return: int, number of active groups
        -----------------------------------------------------------------------------------------"""
        dmat = self.dmat
        active = self.active()

        # update the distances, branch length is stored as 1/2 the distance between the taxa being
        # joined.  This is not quite correct since the length of the children needs to be subtracted
        branch = dmat[row][col] / 2.0
        for i in active:
            # dmat[row][i] = (dmat[row][i] + dmat[col][i]) / 2.0
            dmat[i][row] = (dmat[i][row] + dmat[i][col]) / 2.0
            dmat[row][i] = dmat[i][row]

        node_row = self.tree[active.index(row)]
        node_col = self.tree[active.index(col)]

        # create a new Tree node for the merged taxa
        new_node = Tree()
        new_node.idx = node_row.idx
        new_node.id = f'({node_row.id}:{branch:.3f},{node_col.id}:{branch:.3f})'

        self.tree[row] = new_node
        self.tree.remove(node_col)

        return len(self.tree)


def upgma2(dmat):
    """---------------------------------------------------------------------------------------------
        construct upgma tree

        :param dmat: list of list, symmetric distance matrix
        :return:
    ---------------------------------------------------------------------------------------------"""
    groups = [i for i in range(len(dmat))]
    names = [str(i) for i in range(len(dmat))]

    ngroups = len(groups)
    while ngroups > 1:
        row, col = smallest(dmat, groups)
        branch = dmat[row][col] / 2.0

        names[row] = f'({names[row]}:{branch:.3f},{names[col]}:{branch:.3f})'
        ngroups = mergerow(dmat, groups, row, col)

    return names[groups[0]]


def mergerow(dmat, groups, row, col):
    """---------------------------------------------------------------------------------------------
    merge the row and column with the smallest value

    :param dmat: list of list, symmetric distance matrix
    :param groups: list of int, active indices in dmat
    :param row: index of row (will be kept)
    :param col: index of column (will be merged)
    :return: int, number of active groups
    ---------------------------------------------------------------------------------------------"""
    for i in range(len(dmat[col])):
        dmat[row][i] = (dmat[row][i] + dmat[col][i]) / 2.0

    groups.remove(col)
    return len(groups)


def smallest2(dmat, groups, initval=10):
    """---------------------------------------------------------------------------------------------
    find the smallest value in the current distance matrix, only look at values in groups, other
    indices have been merged

    :param dmat: list of list, symmetric distance matrix
    :param groups: list of int, active indices in dmat
    :return: int, int, indices of the smallest distance
    ---------------------------------------------------------------------------------------------"""
    smallest = initval
    srow = 0
    scol = 0
    for i in range(len(groups)):
        row = groups[i]
        for j in range(i + 1, len(groups)):
            col = groups[j]
            if dmat[row][col] < smallest:
                smallest = dmat[row][col]
                srow = row
                scol = col

    return srow, scol


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # single linkage clusters
    distance_file = 'distance.out'
    threshold = 0.5
    distance, maximum, minimum = read_distance(distance_file)
    cluster, index = connected(distance, threshold)

    # for each cluster, make upgma tree
    for c in range(len(cluster)):
        taxa_n = len(cluster[c])
        if taxa_n > 1:
            tree = Upgma()
            tree.load(distance, cluster[c])
            i = 0
            for taxon in cluster[c]:
                print(f'{i}\t{taxon}')
                i += 1

            print('\n')

            print(tree.dmat_format())
            while taxa_n > 1:
                row, col = tree.smallest()
                taxa_n = tree.mergetaxa(row, col)
                print(tree.active())
                print(tree.dmat_format())

            print(tree.tree[0].id)

    #         print(f'cluster {c}= {cluster[c]}')
    #         dmat = distance_matrix(distance, cluster[c])
    #         print(dmat)
    #         upgma(dmat)

    exit(0)
