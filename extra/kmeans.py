import numpy
import random


class Kmeans():
    """=============================================================================================
    Simple kmeans based on fingerprint
    ============================================================================================="""
    import numpy
    import random

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
        error = 0.0
        error_old = 0.0
        cycle = 0

        centroid = [[0 for _ in range(nfeature)] for _ in range(ngroup)]
        newgroup = [[] for _ in range(ngroup)]

        while True:
            cycle += 1
            maxval = 0

            # calculate centroid positions
            # maxval is the sum of the lengths of the data vectors; used as initial value for mindist
            for g in range(ngroup):
                centroid[g][:] = [0 for g in range(nfeature)]
                for i in group[g]:
                    centroid[g] += self.data[i]
                    # maxval = max(maxval, sum(self.data[i]))
                    maxval += sum(self.data[i])

                if len(group[g]):
                    centroid[g] /= len(group[g])

            maxval *= 2
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

                # reassign point to closest centroid
                try:
                    newgroup[minid].append(point)
                except IndexError:
                    print('oops')

                # add the distance to the current assigned centroid to the error
                error += mindist

                if minid > ngroup:
                    # should never happen
                    # maxval = mindist + 1
                    minid = ngroup - 1

            # copy new groups into group and reset newgroup
            try:
                group = newgroup[:]
            except IndexError:
                print('oops')
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