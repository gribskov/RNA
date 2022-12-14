"""=================================================================================================
Read the distances calculated by fingerprint_distance. Expected format
fpt1    fpt2    Jaccard Bray-Curtis

First cluster in to connected components, then run UPGMA on each cluster

Michael Gribskov     06 May 2022
================================================================================================="""
import sys
from datetime import datetime


def read_distance(filename):
    """---------------------------------------------------------------------------------------------
    Read distances from fingerprint_distance.py.  Distances are stored in a list of dict
    distance = [ {fpt1, fpt2, jaccard, bray_curtis}, ...]. Also counts the number of cases where
    the group of fpt1 matches fpt2 (true positive) and where they don't match (false positive
    
    :param filename: string, name of file with distances
    :return: distance (list of dict), maximum (float), minimum (float), pos(int), neg(int)
    ---------------------------------------------------------------------------------------------"""
    try:
        distance = open(filename, 'r')
    except OSError:
        sys.stderr.write(f'distance_cluster - read_distance: cannot open distance file{filename}{newline}')

    maximum = {'jaccard': 0, 'bray-curtis': 0}
    minimum = {'jaccard': 1000000, 'bray-curtis': 1000000}

    pos = 0
    neg = 0
    ispos = None
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

        f1 = get_group(fpt1)
        f2 = get_group(fpt2)
        ispos = False
        if f1 == f2:
            ispos = True
            pos += 1
        else:
            neg += 1

        distance_list.append({'fpt1':        fpt1,
                              'fpt2':        fpt2,
                              'jaccard':     float(jaccard),
                              'bray-curtis': float(bray_curtis),
                              'ispos':       ispos})

    return distance_list, maximum, minimum, pos, neg


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
            if node is not None:
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

    def similarity_to_distance(self, big):
        """-----------------------------------------------------------------------------------------
        convert a similarity measure (such as Jaccard) to distance by subtracting each value from
        the largest value.

        :param big: int largest value in dmat
        :return:int, int    new min and max of dmat
        -----------------------------------------------------------------------------------------"""
        dmat = self.dmat
        size = len(dmat)
        dmax = dmin = dmat[0][0]
        for i in range(size):
            for j in range(size):
                dmat[i][j] = big - dmat[i][j]
                dmax = max(dmat[i][j], dmax)
                dmin = min(dmat[i][j], dmin)

        return dmin, dmax

    def smallest(self, initval=10):
        """---------------------------------------------------------------------------------------------
        find the smallest value in the current distance matrix, only look at values in groups, other
        indices have been merged

        :return: int, int, indices of the smallest distance in dmat coordinates
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

        # print(smallest, srow, scol)
        return srow, scol

    def mergetaxa(self, row, col):
        """-----------------------------------------------------------------------------------------
        merge the row and column taxa with the smallest value. the row and column are given in dmat
        coordinates

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

        arow = active.index(row)
        acol = active.index(col)
        # node_row = self.tree[active.index(row)]
        # node_col = self.tree[active.index(col)]
        node_row = self.tree[row]
        node_col = self.tree[col]

        # create a new Tree node for the merged taxa
        new_node = Tree()
        new_node.idx = node_row.idx
        new_node.id = f'({node_row.id}:{branch:.3f},{node_col.id}:{branch:.3f})'

        self.tree[row] = new_node
        # print(len(active), row, self.tree[row].id)
        # self.tree.remove(node_col)
        self.tree[col] = None

        return len(active) - 1


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
        row, col = smallest2(dmat, groups)
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


def get_group(name):
    """---------------------------------------------------------------------------------------------
    extract the part of the fingerprint file name up to the first period or underline
    :param name: string fingerprint file name
    :return: string     group name
    ---------------------------------------------------------------------------------------------"""
    breakpoint = len(name)
    for c in ('.', '_'):
        pos = name.find(c)
        if pos > 0:
            breakpoint = min(breakpoint, pos)

    return name[:breakpoint]


def ROC2(roc_file, distance, score, npos, nneg):
    """---------------------------------------------------------------------------------------------
    calculate Receiver Operating Characteristic and area under curve (AUC)
    :param roc_file:    string, filename for ROC output
    :param distance:    list of dict, distance information
    :param npos:        int, number of positive comparisons
    :param nneg:        int, number of negative comparisons
    :return:
    ---------------------------------------------------------------------------------------------"""

    nstep = 1.0 / nneg
    pstep = 1.0 / npos
    pbeg = pend = 0
    nbeg = nend = 0
    dirup = True
    value = 0
    auc = area = 0.0
    for point in sorted(distance, key=lambda x: x[score], reverse=True):
        # print(point)

        if point[score] == value:
            # items with same value cannot be sorted so they are a single step
            if point['ispos']:
                pend += 1
            else:
                nend += 1
            continue

        else:
            # value has changed but you don't need to calculate until both pos
            # and negative values change
            if not dirup and not point['ispos']:
                # this step is negative, following previous negative, just increment neg
                nend += 1
            if dirup and point['ispos']:
                # positive step following previous positives, i.e., straight up
                pend += 1
                pbeg += 1
            else:
                # area of trapezoid
                area = (pbeg + pend) * pstep * (nend - nbeg) * nstep / 2.0
                auc += area
                print(f'score:{point[score]}\t area:{area}\tauc:{auc}')

                if point['ispos']:
                    pend += 1
                    dirup = True
                else:
                    dirup = False
                    nend += 1

                pbeg = pend
                nbeg = nend

                value = point[score]

    area = (pbeg + pend) * pstep * (nend - nbeg) * nstep / 2.0
    auc += area

    print(f'ROC AUC = {auc}')

    return


def ROC(roc_file, distance, score, npos, nneg):
    """---------------------------------------------------------------------------------------------
    calculate Receiver Operating Characteristic and area under curve (AUC)
    :param roc_file:    string, filename for ROC output
    :param distance:    list of dict, distance information
    :param npos:        int, number of positive comparisons
    :param nneg:        int, number of negative comparisons
    :return:
    ---------------------------------------------------------------------------------------------"""
    roc = 'pbegin\tpend\tnbegin\tnend\tarea\tAUC\n'
    auc = area = 0.0

    nstep = 1.0 / nneg
    pstep = 1.0 / npos
    pbeg = pend = 0
    nbeg = nend = 0
    value = 1.0

    # print(f'{npos}\t{nneg}')
    for point in sorted(distance, key=lambda x: x[score], reverse=True):
        # print(f'{pbeg}\t{pend}\t{nbeg}\t{nend}\t{point["jaccard"]}\t{point["ispos"]}')

        # detect blocks where items have the same value, these cannot be sorted
        # so they are a single step
        if point[score] == value:
            if point['ispos']:
                pend += pstep
            else:
                nend += nstep
            continue

        else:
            # we get here if we have completed a block of equal valuesif both p and n ranges
            # calculate trapezoidal area
            area = (pbeg + pend) * (nend - nbeg) / 2.0
            auc += area
            roc += f'{pbeg:.3g}\t{pend:.3g}\t{nbeg:.3g}\t{nend:.3g}\t{point[score]:.3g}\t{area:.3g}\t{auc:.3g}\n'
            nbeg = nend
            pbeg = pend
            if point['ispos']:
                pend += pstep
            else:
                nend += nstep

            value = point[score]

    area = (pbeg + pend) * (nend - nbeg) / 2.0
    auc += area
    print(f'{pbeg:.3g}\t{pend:.3g}\t{nbeg:.3g}\t{nend:.3g}\t{point[score]:.3g}\t area:{area:.3g}\tauc:{auc:.3g}')

    return roc, auc

def make_trees(opt, distance, cluster):
    """---------------------------------------------------------------------------------------------

    :param opt:
    :param distance:
    :param cluster:
    :return:
    ---------------------------------------------------------------------------------------------"""


def process_command_line():
    """---------------------------------------------------------------------------------------------
    read command line options and return as a Namespace object. The Namespace object behaves much as
    a dictionary, and can be converted to a dictionary using the vars() method

    :return: Namespace object
    ---------------------------------------------------------------------------------------------"""
    import argparse

    cl = argparse.ArgumentParser(
        description='Cluster and calculate ROC from distances',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=140, max_help_position=80)
        )
    cl.add_argument('-t', '--type',
                    help='distance type: jaccard or bray-curtis (default=%(default)s)', type=str,
                    choices=['jaccard', 'bray-curtis'], nargs='*',
                    default=['jaccard'])
    # action = 'extend',
    cl.add_argument('-d', '--distance',
                    help='distance file from fingerprint_distance (default=%(default)s)',
                    default='fingerprint.distance')
    cl.add_argument('-m', '--mindist',
                    help='Minimum distance for clustering (default=%(default)s)', type=float,
                    default=0.0)


    rocoptions = cl.add_argument_group('ROC options')
    rocoptions.add_argument('-r', '--roc',
                    help='do not calculate ROC and AUC', action='store_false')
    rocoptions.add_argument('-p', '--plotable', action='store_true',
                    help='save ROC histogram for plotting')

    treeoptions = cl.add_argument_group('Tree options')
    treeoptions.add_argument('-c', '--condensed', action='store_false',
                    help='include Newick formatted tree (default=%(default)s)',
                    default=False)
    treeoptions.add_argument('-i', '--indented', action='store_false',
                    help='include indented Newick formatted tree (default=%(default)s)',
                    default=False)

    # TODO cluster prefixes for tree output (tree output)
    # use False for none? or add tree option prefix, stdout

    args = cl.parse_args()

    return args


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    newline = '\n'

    opt = process_command_line()
    sys.stderr.write(f'distance_cluster: Cluster and calculate ROC from distance{newline}')
    sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{newline}{newline}")
    sys.stderr.write(f"distance file: {opt.distance}{newline}")
    sys.stderr.write(f"clustering threshold: {opt.mindist}{newline}")

    # single linkage clusters
    distance, maximum, minimum, pos, neg = read_distance(opt.distance)

    if opt.roc:
        roc_file = 'roc.out'
        roc, auc = ROC(roc_file, distance, 'jaccard', pos, neg)
        sys.stderr.write(f'{newline}ROC AUC = {auc:.4g}{newline}')

    cluster, index = connected(distance, opt.mindist)

    # for each cluster, make upgma tree
    if opt.condensed or opt.indented
    for c in range(len(cluster)):
        taxa_n = len(cluster[c])
        sys.stdout.write(f'# Cluster_{c}: {taxa_n} fingerprints')

        if taxa_n > 0:
            tree = Upgma()

            tree.load(distance, cluster[c])
            if tree.dtype == 'jaccard':
                minimum['jaccard'], maximum['jaccard'] = tree.similarity_to_distance(maximum['jaccard'])

            i = 0
            print(f'cluster: {c} input taxa')
            for taxon in cluster[c]:
                print(f'{i}\t{taxon}')
                i += 1

            while taxa_n > 1:
                row, col = tree.smallest()
                taxa_n = tree.mergetaxa(row, col)

            # final tree is always taxon 0
            print(f'\nfinal tree')
            tree.tree[0].id += ';'
            print(tree.tree[0].id)

            # indented Newick pseudo tree
            print(f'\nindented tree')
            pos = 0
            tree = tree.tree[0].id
            indent = 0
            start = 0
            while pos < len(tree):
                if tree[pos] == '(':
                    print(f'{" " * indent}(')
                    indent += 4
                    start = pos + 1
                elif tree[pos] == ')':
                    print(f'{" " * indent}{tree[start:pos]}')
                    indent -= 4
                    print(f'{" " * indent})', end='')
                    start = pos + 1

                elif tree[pos] == ',':
                    if tree[start] == ':':
                        print(f'{tree[start:pos]},')
                    else:
                        print(f'{" " * indent}{tree[start:pos]},')
                    start = pos + 1

                pos += 1

    exit(0)
