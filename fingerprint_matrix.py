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

    def cluster(self):

        # get cluster centroids
        nfeature = len(self.data[0])
        for group in self.group:
            self.centroid.append( [0 for _ in range(nfeature)])
            for i in group:
                self.centroid[-1] += self.data[i]
            self.centroid[-1] /= len(group)

        for g in range(len(self.group)):
            for point in self.group[g]:
                dist = sum(abs(self.data[point] - self.centroid[g]))

        return

# main program
if __name__ == '__main__':
    selection = 'data/fpt/*'
    fmat = FingerprintMatrix()
    fmat.read_files(selection)
    # fmat.select_min_max(2, 2)
    # fmat.select_min_max(0, 7, setval=True, recalculate=True)
    # fmat.index2matrix()
    k = Kmeans(2)
    ndata = k.assign_data(fmat.fpt)
    print(f'{ndata} points read for k={k.k}')

    k.cluster()

    exit(0)
