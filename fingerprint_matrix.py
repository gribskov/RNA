"""=================================================================================================
Convert a set of fingerprint files into a binary form for rapid reading and manipulation
================================================================================================="""
from fingerprint import FingerprintMatrix
import numpy
from math import floor, ceil

import random


class Kmeans():
    """=============================================================================================
    Simple kmeans based on fingerprint
    ============================================================================================="""

    def __init__(self, k=2):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.k = k
        self.group = None
        self.data = None
        self.centroid = []
        self.error = 0

    def assign_data(self, mat):
        """-----------------------------------------------------------------------------------------


        :param mat:list of list     feature vector data
        :return:
        -----------------------------------------------------------------------------------------"""
        # distribute indices to k groups
        group = self.group
        k = self.k
        group = [[] for _ in range(self.k)]
        i = 0
        for i in range(len(mat)):
            g = i % k
            group[g].append(i)

        self.group = group
        for i in range(len(mat)):
            mat[i] = mat[i].astype(float)

        self.data = mat

        return len(self.data)

    def assign_data_random(self, mat):
        """-----------------------------------------------------------------------------------------
        Randomly choose k datapoint to be initial centroids

        :param mat:list of list     feature vector data
        :return:
        -----------------------------------------------------------------------------------------"""
        # distribute indices to k groups
        k = self.k

        index = [i for i in range(len(mat))]
        random.shuffle(index)

        self.group = []
        for g in range(k):
            self.group.append([index[g]])

        for i in range(len(mat)):
            # convert true false in data vectors to float
            mat[i] = mat[i].astype(float)

        self.data = mat

        return len(self.data)

    def cluster(self, show_cycle=False, delta=1.0e-4):
        """-----------------------------------------------------------------------------------------
        perform one round of clustering, cycle until the change in error is < delta

        :param show_cycle: bool     show clusters at each cycle
        :param delta: float         error cutoff for termination
        :return: bool               True for now
        -----------------------------------------------------------------------------------------"""
        numpy.seterr(all='raise')

        # get cluster centroids
        group = self.group
        centroid = self.centroid
        nfeature = len(self.data[0])
        npoint = len(self.data)
        ngroup = len(self.group)
        maxval = 0
        error = 0.0
        error_old = 0.0
        cycle = 0

        centroid = [[0 for _ in range(nfeature)] for _ in range(ngroup)]
        newgroup = [[] for _ in range(ngroup)]

        while True:
            cycle += 1

            # calculate centroid positions
            for g in range(ngroup):
                centroid[g][:] = [0 for g in range(nfeature)]
                for i in group[g]:
                    centroid[g] += self.data[i]
                    # maxval = max(maxval, sum(self.data[i]))
                    maxval += sum(self.data[i])

                if len(group[g]):
                    centroid[g] /= len(group[g])

            error = 0
            for point in range(npoint):
                # for each point find the distance to each centroid and reassign
                mindist = maxval
                minid = ngroup + 1
                for g in range(ngroup):
                    # each current group centroid
                    dist = sum(abs(self.data[point] - centroid[g]))
                    # dist = sum(abs(self.data[point] - centroid[g])) / max(sum(self.data[point]), 1)
                    # dist = sum(abs(self.data[point] - centroid[g])) * max(sum(self.data[point]), 1)
                    if dist <= mindist:
                        mindist = dist
                        minid = g

                # reassing point to closest centroid
                newgroup[minid].append(point)

                # add the distance to the current assigned centroid to the error
                error += mindist

                if minid > ngroup:
                    maxval = mindist + 1
                    minid = ngroup - 1

            # copy new groups into group and reset newgroup
            group = newgroup[:]
            for g in range(ngroup):
                newgroup[g] = []

            # current error
            error /= npoint
            error_relative = error
            if error > 0:
                error_relative = abs(error_old - error) / error
            error_old = error

            print(f'\tcycle={cycle}\t error={error:.1f} rel.error={error_relative:.3g}')
            if show_cycle:
                for g in range(ngroup):
                    print(f'\t{g}\t{group[g]}')

            if error_relative < delta:
                break

        self.centroid = centroid
        self.group = group

        return True


###################################################################################################
# main program
###################################################################################################
if __name__ == '__main__':
    n_repeat = 150
    k = Kmeans(30)
    neighborhood_cutoff = (0.4, 0.6)
    show_cycle = False
    show_init = False

    selection = 'data/fpt/*.out'
    fmat = FingerprintMatrix()
    fmat.read_files(selection)
    fmat.select_min_max(10, 45, False, recalculate=True)


    together = [[0 for _ in range(len(fmat.fpt))] for _ in range(len(fmat.fpt))]

    for repeat in range(n_repeat):
        # cluster n_repeat times
        # TODO save the lowest error clustering
        ndata = k.assign_data_random(fmat.fpt)
        if show_init:
            print(f'{ndata} points from {selection} for k={k.k}')
            print(f'initial groups')
            for g in range(len(k.group)):
                print(f'group {g}: ', end='')
                for i in k.group[g]:
                    print(f'\t{i}', end='')
                print()

        print(f'repeat {repeat}')
        k.cluster()

        g = 0
        fptname = list(fmat.fpt_id.keys())
        for group in k.group:
            g += 1
            if show_cycle:
                # print the clusters at each cycle
                print(f'group {g}')
                for f in group:
                    print(f'\t{fptname[f]}')

            pairs = ((a, b) for idx, a in enumerate(group) for b in group[idx + 1:])
            for a, b in pairs:
                together[a][b] += 1
                together[b][a] += 1

        # final clusters
        # g=0
        # for group in k.group:
        #     g += 1
        #     print(f'group {g}')
        #     for f in group:
        #         print(f'\t{fptname[f]}')

    # neighborhood clustering: group the points that co-occur the most often
    # for minval in range(floor(n_repeat * neighborhood_cutoff), n_repeat):
    for minval in range(floor(n_repeat * neighborhood_cutoff[0]), ceil(n_repeat * neighborhood_cutoff[1])):
        print(f'\n{"#" * 80}\n minval={minval}\n{"#" * 80}')
        cluster = []
        n = 0
        avail = [i for i in range(len(fmat.fpt_id))]
        while avail:
            c = avail.pop()
            cluster.append([c])
            todo = [c]
            while todo:
                p0 = todo.pop()
                for p1 in avail:
                    if together[p0][p1] > minval:
                        cluster[n].append(p1)
                        todo.append(p1)
                        avail.remove(p1)

            n += 1

        c = 0
        for neighborhood in cluster:
            print(f'cluster {c}')
            c += 1
            for f in neighborhood:
                print(f'\t{fptname[f]}')

exit(0)
