"""=================================================================================================
Convert a set of fingerprint files into a binary form for rapid reading and manipulation
================================================================================================="""
from fingerprint import FingerprintMatrix
import numpy


# import random

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

    def cluster(self, delta=0.1):

        # get cluster centroids
        group = self.group
        centroid = self.centroid
        nfeature = len(self.data[0])
        npoint = len(self.data)
        ngroup = len(self.group)
        maxval = 0
        error = 0
        error_old = 0
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
                    maxval = max(maxval, sum(self.data[i]))

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
                    if dist < mindist:
                        mindist = dist
                        minid = g
                    if point in group[g]:
                        # add the distance to the current assigned centroid to
                        # the error
                        error += dist

                newgroup[minid].append(point)

            # current error
            error /= npoint
            print(f'cycle={cycle}\t error={error:.1f}')
            for g in range(ngroup):
                print(f'\t{g}\t{group[g]}')

            group = newgroup[:]
            for g in range(ngroup):
                newgroup[g] = []

            if (abs(error_old - error) / error < delta):
                break

            error_old = error

        self.centroid = centroid
        self.group = group

        return


# main program
if __name__ == '__main__':
    selection = 'data/fpt/oldfpt/*.xpt'
    fmat = FingerprintMatrix()
    fmat.read_files(selection)
    # fmat.select_min_max(2, 2)
    # fmat.select_min_max(0, 7, setval=True, recalculate=True)
    # fmat.index2matrix()
    k = Kmeans(5)
    ndata = k.assign_data(fmat.fpt)
    print(f'{ndata} points read for k={k.k}')

    k.cluster()
    g = 0
    fptname = list(fmat.fpt_id.keys())
    for group in k.group:
        g += 1
        print(f'group {g}')
        for f in group:
            print(f'\t{fptname[f]}')

    exit(0)
