import sys
import glob
import os
import json
import datetime
import yaml
import numpy
import pickle

# note install PyYAML

# for use in fstrings
newline = '\n'


class Fingerprint(dict):
    """=============================================================================================
    A fingerprint is a dict tabulating the spectrum of fixed size motifs in a structure.  The keys
    of the dictionary are the minimum dfs codes of the motifs.  These correspond to the keys of the
    motifs created by make_xios_db.py.

    Sorting the results

    # sort by count
    for motif,count in sorted(fingerprint.motif.items(),key=lambda x:x[1],reverse=True):

    # sort alphabetically
    for motif,count in sorted(fingerprint.motif.items(),key=lambda x:x[0]):

    10 July 2019     Michael Gribskov
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        The main data structure is the dictionary of motifs.  Normally the keys are the minimum DFS
        and the values are the counts
        -----------------------------------------------------------------------------------------"""
        super().__init__(self)
        self.information = {'Date': ''
                            }
        self.motif = {}
        self.count = 0  # sum of counts of all motifs

        self.setdate()

    @property
    def n(self):
        """-----------------------------------------------------------------------------------------
        makes a function for length that looks like an attribute (but read only)
        e.g., print(fp.n)

        :return: int, number of motifs in motif dictionary
        -----------------------------------------------------------------------------------------"""
        return len(self.motif)

    def add(self, string, n=1):
        """-----------------------------------------------------------------------------------------
        Adds a motif to the fingerprint, creating a new entry if needed.  Agnostic as to the
        encoding of the motif string

        :param string: str, encoded motifs string
        :param n: int, number of observations of this motif
        :return: int, number of motifs in fingerprint
        -----------------------------------------------------------------------------------------"""
        self.count += n
        if string in self.motif:
            self.motif[string] += n
        else:
            self.motif[string] = n
            self.count += n

        return self.count

    def get(self, key):
        """-----------------------------------------------------------------------------------------
        Lookup and return a motif by name (key)

        :param key: string, a key in self.motif
        :return:
        -----------------------------------------------------------------------------------------"""
        if key in self.motif:
            return self.motif[key]
        else:
            return None

    def mincount(self):
        """-----------------------------------------------------------------------------------------
        Return the count of the motif with the smallest count

        :return: int
        -----------------------------------------------------------------------------------------"""
        return min(self.motif.values(), default=0)

    def minkey(self):
        """-----------------------------------------------------------------------------------------
        Return the motif with the lowest count

        :return: str
        -----------------------------------------------------------------------------------------"""
        m = self.motif
        m = min(m.items(), default=0, key=lambda x: x[1])
        # print('smallest', m)

        return m[0]

    def toJSON(self):
        """-----------------------------------------------------------------------------------------
        Convert fingerprint to JSON string

        :return: str
        -----------------------------------------------------------------------------------------"""
        fields = ['information', 'total', 'nmotif', 'motif']

        root = []
        for i in range(len(fields)):
            # look up the values of each field and add to dispatch list
            root.append(getattr(self, fields[i]))

        return json.dumps({fields[i]: root[i] for i in range(len(fields))}, indent=4)

    def toYAML(self, sort='len'):
        """-----------------------------------------------------------------------------------------
        Convert fingerprint to JSON string

        :return: str
        -----------------------------------------------------------------------------------------"""
        fields = ['information', 'total', 'nmotif', 'motif']

        m = self.motif
        if sort == 'len':
            m = {}
            for k in sorted(self.motif, key=lambda x: self.motif[x], reverse=True):
                m[k] = self.motif[k]

        root = [{'fingerprint': [{'information': self.information},
                                 {'total': self.count},
                                 {'nmotif': self.n},
                                 {'motif': m}]
                 }]

        return yaml.dump(root, indent=2, default_flow_style=False, sort_keys=False)
        # return yaml.dump(root, indent=2, default_flow_style=False)

    def writeYAML(self, file):
        """-----------------------------------------------------------------------------------------
        Write the fingerprint to a file in YAML format

        :param file: fp/str, either a string or an open file
        :return: True
        -----------------------------------------------------------------------------------------"""
        if isinstance(file, str):
            # file argument is string, try to open
            try:
                fp = open(file, 'w')
            except OSError:
                sys.stderr.write('fingerprint.writeYAML - error opening file ({})'.format(file))
                exit(1)
        else:
            # file is not str, assume it is a file pointer
            fp = file

        fp.write(self.toYAML())
        # fp.close()    if fp is stdout, not good

        return True

    def readYAML(self, file):
        """-----------------------------------------------------------------------------------------
        read the fingerprint from a file in YAML format

        :param file: fp/str, either a string or an open file
        :return: int, number of motfs
        -----------------------------------------------------------------------------------------"""
        if isinstance(file, str):
            # file argument is string, try to open
            try:
                fp = open(file, 'r')
            except OSError:
                sys.stderr.write('fingerprint.readYAML - error opening file ({})\n'.format(file))
                exit(1)
        else:
            # file is not str, assume it is a file pointer
            fp = file

        f = yaml.load(fp, Loader=yaml.FullLoader)
        if f == None:
            sys.stderr.write('No fingerprint found in {}\n'.format(file))
            self.information = {'File': file}
            self.count = 0
            self.motif = {}
        else:
            # fields = ['information', 'total', 'nmotif', 'motif']
            root = f[0]['fingerprint']
            self.information = root[0]['information']
            self.count = root[1]['total']
            self.motif = root[3]['motif']

        return self.n

    def setdate(self):
        """-----------------------------------------------------------------------------------------
        set date in information
        
        :return: datetime
        -----------------------------------------------------------------------------------------"""
        daytime = datetime.datetime.now()
        self.information['Date'] = daytime.strftime('%Y-%m-%d %H:%M:%S')

        return daytime

    def add_parents(self, motifdb):
        """-----------------------------------------------------------------------------------------
        Look up the parents in the provided motif database and add them to the fingerprint.  The
        count of each parent is incremented for each of its children.  for nested children this may
        not make perfect sense

        :param motifdb: MotifDB object with motifs and parents
        :return: int, total number of motifs
        -----------------------------------------------------------------------------------------"""
        children = list(self.motif.keys())
        # print(f'children:{self.n}')
        i = 0
        for child in children:
            i += 1
            try:
                for parent in motifdb.parent[child]:
                    # print('child {}\t\tparent {} {}'.format(child, parent, self.motif[child]))
                    self.add(parent, n=self.motif[child])
            except KeyError:
                # shouldn't see this, of course
                print(f'{i}\tadd_parents error:{child}')

        return self.n


class FingerprintSet(list):
    """=============================================================================================
    A collection of fingerprints and similarity/distance functions

    TODO: instead of having a long list of similarity/distance function, add a dispatcher that
    handles both similarity and distance for all implemented types (and also sets precision)
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Fingerprint set is basically just a list with extra methods, since fingerprint inherits
        from list, self is a list.

        in addition to converting between numeric and string indices of motifs, i2motif give a list
        of the names of all motifs
        -----------------------------------------------------------------------------------------"""
        super().__init__(self)
        self.motif2i = {}
        self.i2motif = []

    def fill(self):
        """-----------------------------------------------------------------------------------------
        compare the motif lists for all the fingerprint in the set and fill in any missing  motifs
        with zeroes.  After running fill, all fingerprints will have the same motif lists

        :return: int, number of motifs
        -----------------------------------------------------------------------------------------"""
        # first pass to get a list of all motifs
        motiflist = []
        if len(self) == 0:
            return 0

        first = True
        for fingerprint in self:
            if first:
                motiflist = self[0].keys()
                first = False
                continue

            for motif in fingerprint.motif:
                if not motif in motiflist:
                    motiflist.append(motif)

        # second pass to fill in missing values
        for fingerprint in self:
            for motif in motiflist:
                if not motif in fingerprint.motif:
                    fingerprint.motif[motif] = 0

        return len(motiflist)

    def jaccard_binary(self):
        """-----------------------------------------------------------------------------------------
        calculate jaccard similarity using the binary motif matrix
        TODO is this faster in numpy?
        :return:
        -----------------------------------------------------------------------------------------"""
        if not self.motif_matrix:
            # if motif index does not exist, use all motifs, for selected motifs, index is made by
            # index_all_motifs()
            self.index_all_motifs()

        mm = self.motif_matrix
        j_sim = []
        for i in range(0, len(self)):
            for j in range(i + 1, len(self)):
                intersect = 0
                union = 0
                vi = mm[i]
                vj = mm[j]
                for k in range(0, len(self.i2motif)):
                    intersect += vi[k] and vj[k]
                    union += vi[k] or vj[k]

                try:
                    jaccard = intersect / union
                except ZeroDivisionError:
                    jaccard = 0

                j_sim.append([i, j, jaccard])

        return j_sim

    def jaccard_scale(self):
        """-----------------------------------------------------------------------------------------
        calculate jaccard similarity using the binary motif matrix
        TODO is this faster in numpy?
        :return:
        -----------------------------------------------------------------------------------------"""
        if not self.motif_matrix:
            # if motif index does not exist, use all motifs, for selected motifs, index is made by
            # index_all_motifs()
            self.index_all_motifs()

        mm = self.motif_matrix
        j_sim = []
        for i in range(0, len(self)):
            for j in range(i + 1, len(self)):
                intersect = 0
                union = 0
                vi = mm[i]
                vj = mm[j]
                for k in range(0, len(self.i2motif)):
                    intersect += vi[k] and vj[k]
                    # union += vi[k] or vj[k]
                minfpt = min(sum(vi), sum(vj))

                try:
                    jaccard = intersect / minfpt
                except ZeroDivisionError:
                    jaccard = 0

                j_sim.append([i, j, jaccard])

        return j_sim

    def jaccard_sim(self, idx=[]):
        """-----------------------------------------------------------------------------------------
        Calculate pairwise Jaccard similarity between the fingerprints indicated by idx.  [0,1],
        e.g., means just the first two fingerprints in the set.

        range: [0,1]

        :param idx: list, indices of members of FingerprintSet to compare,
        :return: list of float, similarity values
        -----------------------------------------------------------------------------------------"""
        nfp = len(self)
        if len(idx) == 0:
            idx = [i for i in range(nfp)]
        elif max(idx) > len(self) - 1:
            sys.stderr.write(f'FingerprintSet:jaccard - selected indices exceed number of motifs\n')
            exit(2)

        jaccard = []
        for i in range(len(idx)):
            m_i = self[idx[i]].motif

            for j in range(i + 1, len(idx)):
                m_j = self[idx[j]].motif
                intersect = []
                union = []

                for motif in m_i:
                    if motif in m_j:
                        intersect.append(motif)
                    union.append(motif)

                for motif in m_j:
                    if motif not in union:
                        union.append(motif)

                try:
                    jaccard.append([idx[i], idx[j], len(intersect) / len(union)])
                except ZeroDivisionError:
                    sys.stderr.write(
                        f'FingerprintSet:jaccard - no motifs in fingerprints in {idx[i]} and {idx[j]}\n')

        return jaccard

    def bray_curtis_dis(self, idx=[]):
        """-----------------------------------------------------------------------------------------
        Calculate pairwise Bray-Curtis dissimilarity between the fingerprints indicated by idx.
        [0,1], e.g., means just the first two fingerprints in the set.

        range: [0,1]

        :param idx: list, indices of members of FingerprintSet to compare,
        :return: list of float, dissimilarity values
        -----------------------------------------------------------------------------------------"""
        nfp = len(self)
        if len(idx) == 0:
            idx = [i for i in range(nfp)]
        elif max(idx) > len(self) - 1:
            sys.stderr.write(
                f'FingerprintSet:bray_curtis_dis - selected indices exceed number of motifs\n')
            exit(2)

        bc = []
        for i in range(len(idx)):
            m_i = self[idx[i]].motif

            for j in range(i + 1, len(idx)):
                m_j = self[idx[j]].motif

                intersect = 0
                union = 0
                for motif in m_i:
                    # print(f'bc {i} x {j}')
                    if motif in m_j:
                        intersect += min(m_i[motif], m_j[motif])
                    union += m_i[motif]

                for motif in m_j:
                    union += m_j[motif]

                try:
                    bc.append([idx[i], idx[j], 2 * intersect / union])
                except ZeroDivisionError:
                    sys.stderr.write(
                        f'FingerprintSet:bray_curtis_dis - no motifs in fingerprints in {idx[i]} and'
                        f' {idx[j]}\n')

        return bc

    def bray_curtis_binary(self):
        """-----------------------------------------------------------------------------------------
        Calculate pairwise Bray-Curtis dissimilarity between the fingerprints based on the binary
        motif matrix. This is not precisely Bray-Curtis which would use the counts of the motifs

            BC = 1 - 2 * common_motifs / (n_motifs_i + n_motifs_j)

        :return: list of float, dissimilarity values
        -----------------------------------------------------------------------------------------"""
        nfp = len(self)
        if not self.motif_matrix:
            # if motif index does not exist, use all motifs, for selected motifs, index is made by
            # index_all_motifs()
            self.index_all_motifs()

        mm = self.motif_matrix
        braycurtis = []
        for i in range(nfp):
            vi = mm[i]
            i_n = sum(vi)
            if i_n == 0:
                sys.stderr.write(
                    f'FingerprintSet:bray_curtis_binary - no motifs in fingerprint {i} ')
                sys.stderr.write(f'({self.i2motif[i]})\n')
                for j in range(i + 1, nfp):
                    braycurtis.append([i, j, 1.0])

                # go on to next i
                break

            for j in range(i + 1, nfp):
                intersect = []
                union = []
                vj = mm[j]
                j_n = sum(vj)
                intersect = 0
                union = 0
                for k in range(0, len(self.i2motif)):
                    intersect += vi[k] and vj[k]
                    union += vi[k] or vj[k]

                try:
                    bc = 1.0 - 2.0 * intersect / (i_n + j_n)
                except ZeroDivisionError:
                    bc = 1

                braycurtis.append([i, j, bc])

        return braycurtis

    def select(self, filename=None, selected=None):
        """-----------------------------------------------------------------------------------------
        Make a list of all motifs to be used, and translation to and from numeric indices
        (motif2i, i2motif). If a selected list is not provided, use all motifs.

        :param filename: string     optional filename for an external file with selected motifs
        :param selected: list       motif names, e.g. 0o1.1o2.2o3.3o4.
        :return: int                number of motifs in original set of fingerprints
        -----------------------------------------------------------------------------------------"""
        # check to see if a file of selected motifs is provided
        if os.access(filename, os.R_OK):
            # the file can be opened, so read it
            motifs = open(filename, 'r')
            for line in motifs:
                if line.startswith('#'):
                    # selected motifs file can have comments
                    continue
                mname, count = line.rstrip().split()
                self.i2motif.append(mname)

            motifs.close()
            motif_n = len(self.i2motif)

        else:
            # selected motifs file couldn't be opened (or is absent)
            if filename:
                sys.stderr.write(f'\nmotif file {filename} is not accessible\n')
                sys.stderr.write(f'All motifs will be used\n')
            motif_n = self.index_all_motifs()

        # if a list of selected motifs is provided, shorten the selection list to just the ones
        # provided. motifs that are not present in the motif list above are skipped
        if selected:
            new_list = []
            for motif in selected:
                if motif in self.i2motif:
                    new_list.append(motif)

            self.i2motif = new_list
            # do not update motif_n to the new list

        # set up translation between motif id and numeric index
        if selected:
            # selected motifs
            for i in range(motif_n):
                self.motif2i[self.i2motif[i]] = i

        return motif_n

    def binary_matrix(self):
        """-----------------------------------------------------------------------------------------
        For selected motifs (see select() make a binary matrix indicating presence of absence of
        each motif
        -----------------------------------------------------------------------------------------"""
        motif2i = self.motif2i
        i2motif = self.i2motif
        motif_matrix = []
        for fpt in self:
            motif_vector = [0 for _ in range(len(self.i2motif))]
            motif_matrix.append(motif_vector)
            for motif in fpt.motif:
                if motif in i2motif:
                    motif_vector[i2motif.index(motif)] = 1

        self.motif_matrix = motif_matrix

    def index_all_motifs(self):
        """-----------------------------------------------------------------------------------------
        make an index of all the motifs found in the fingerprintSet.
        motif2i = dict keys=motif IDs, values = numeric index
        i2motif = list, index is numeric index of motif, value is ID string

        :return: int number of unique motifs
        -----------------------------------------------------------------------------------------"""
        motif2i = self.motif2i
        i2motif = self.i2motif

        for fpt in self:
            for motif in fpt.motif:
                if motif not in motif2i:
                    motif2i[motif] = len(i2motif)
                    i2motif.append(motif)

        return len(i2motif)


class FingerprintMatrix:
    """=========================================================================================
    Reading and writing large sets of fingerprints is time consuming. this class creates a
    presence/absence matrix for fingerprints from a group of fingerprint files
    ========================================================================================="""

    def __init__(self):
        """-------------------------------------------------------------------------------------
        motifs has fields id, index, count, selected
        fpt is a dict of list, id, fpt_vector where fpt_vector is the presence/absence vector

        -------------------------------------------------------------------------------------"""
        self.motifs = {}
        self.fpt_id = {}  # index of fingerprints
        self.fpt = []  # list of fingerprints

    def read_files(self, select_str):
        """-------------------------------------------------------------------------------------
        Read selected files and convert to a matrix of true/false indicating presence absence of
        each motif in the set. rows=fingerprints, columns=motifs

        :param select_str:str       file path for fingerprint files, e.g. data/*.fpt
        :return: bool               True
        -------------------------------------------------------------------------------------"""
        motif_n = len(self.motifs)

        fpt_list = glob.glob(select_str)
        #sys.stderr.write(f'{select_str} =>\n\t{fpt_list}\n')
        sys.stderr.write('\n')
        motif_read = 0
        for fpt_file in fpt_list:
            motif_read += 1
            sys.stderr.write(f'\t{motif_read:-3d}  reading {fpt_file} ...\n')
            f = Fingerprint()
            f.readYAML(fpt_file)
            self.fpt_id[fpt_file] = len(self.fpt)
            self.fpt.append([])

            for motif in f.motif:
                if motif in self.motifs:
                    # known motif
                    self.motifs[motif]['count'] += 1

                else:
                    # new motif
                    self.motifs[motif] = {'count': 1, 'index': motif_n, 'selected': True}
                    motif_n += 1

                index = self.motifs[motif]['index']
                self.fpt[-1].append(index)

        self.index2matrix()
        return True

    def motifs_selected(self):
        """-----------------------------------------------------------------------------------------
        Return matrix of selected motif names corresponding to the current fpt matrix

        :return: list of str    motif names
        -----------------------------------------------------------------------------------------"""
        return  [ id for id in self.motifs if self.motifs[id]['selected'] ]


    def write(self, filename):
        """-----------------------------------------------------------------------------------------
        write the matrix to a tab delimited files with 1/0 indicating present/absent
        :param filename:
        :return:
        -----------------------------------------------------------------------------------------"""
        out = None
        try:
            out = open(filename, 'w')
        except OSError:
            sys.stderr.write(f'FingerprintMatrix::write - unable to open file for writing ({filename}\n')
            exit(1)

        # column labels
        for id in self.fpt_id:
            out.write(f'\t{id}')
        out.write('\n')

        # rows are motifs presence/absence = 1/0
        m = 0
        for motif in self.motifs:
            out.write(f'{motif}')
            for f in self.fpt:
                # one row of presence/absence across all fingerprints
                if f[m]:
                    out.write('\t1')
                else:
                    out.write('\t0')

            out.write('\n')
            m += 1

        return True



    def index2matrix(self):
        """-------------------------------------------------------------------------------------
        Convert fingerprints in the form of lists of indices to binary, 1=presence, 2=absence
        indices are in self.motifs
        This method works of self.fpt as a list of indices

        :return:
        -------------------------------------------------------------------------------------"""
        motifs = self.motifs
        # for each fingerprint, create a vector of True/False indicating presence/absence of motif
        # make a vector of selected indices
        indices = [motifs[m]['index'] for m in motifs if motifs[m]['selected']]

        for f in range(len(self.fpt)):
            # make True/False vector for each fingerprint
            iterable = [True if i in self.fpt[f] else False for i in indices]
            # iterable = [True if self.fpt[f][i] else False for i in indices]
            binary = numpy.fromiter(iterable, bool)
            # print(sum(binary))

            self.fpt[f] = binary

        return True

    def update_matrix(self):
        """-------------------------------------------------------------------------------------
        Convert fingerprints in the form of lists of indices to binary, 1=presence, 2=absence
        indices are in self.motifs
        This function works on the binary matrix form of self.fpt

        :return:
        -------------------------------------------------------------------------------------"""
        motifs = self.motifs
        # for each fingerprint, create a vector of True/False indicating presence/absence of motif
        # make a vector of selected indices
        indices = [motifs[m]['index'] for m in motifs if motifs[m]['selected']]

        for f in range(len(self.fpt)):
            # make True/False vector for each fingerprint
            # iterable = [True if i in self.fpt[f] else False for i in indices]
            iterable = [True if self.fpt[f][i] else False for i in indices]
            binary = numpy.fromiter(iterable, bool)
            # print(sum(binary))

            self.fpt[f] = binary

        return True

    def select_min_max(self, minval=0, maxval=None, setval=False, recalculate=False):
        """-----------------------------------------------------------------------------------------
        select motifs that have count >= minval and count <= maxval
        will not change selection status of motifs between minval and maxval
        rebuild the binary matrix if recalculate is True

        :param minval:int           motifs with count<minval are set to setval
        :param maxval:int           motifs with count>maxval are set to setval
        :param setval:bool          value to set motifs < minval or > maxval to
        :param recalculate:bool     if True recalculate the motif matrix
        :return:int                 number of selected motifs
        -----------------------------------------------------------------------------------------"""
        motifs = self.motifs
        motif_n = len(self.fpt)
        if maxval == None:
            maxval = len(self.fpt)

        n_set = 0
        for m in motifs:
            if motifs[m]['count'] < minval or motifs[m]['count'] > maxval:
                motifs[m]['selected'] = setval
                n_set += 1

        if recalculate:
            # TODO index2matrix won't work on the binary matrix,
            # self.index2matrix()
            self.update_matrix()


        return len(self.fpt[0])

    def pickle(self, outfilename):
        """-----------------------------------------------------------------------------------------
        pickle the object and store in outfilename
        :param outfilename:
        :return:
        -----------------------------------------------------------------------------------------"""
        picklefile = open(outfilename, 'wb')
        pickle.dump(self, picklefile)
        return True

    @classmethod
    def unpickle(cls, infilename):
        """-----------------------------------------------------------------------------------------
        reload the pickled object from infilename
        :param infilename:
        :return:
        -----------------------------------------------------------------------------------------"""
        picklefile = open(infilename, 'rb' )
        fmat = pickle.load(picklefile)

        return fmat

    # ##################################################################################################
    # Testing
    # ##################################################################################################
    if __name__ == '__main__':
        finger1 = Fingerprint()
        finger1.add('dummy')
        finger1.add('dummier')
        for i in range(5):
            finger1.add('dummiest')

        key = 'dummier'
        print(f'{key} count={finger1.get(key)}')

        print(f'finger: {finger1}, n: {finger1.n}')
        print(finger1.toYAML())

        finger2 = Fingerprint()
        finger2.add('dummier', n=5)
        finger2.add('x')
        finger2.add('y')

        finger3 = Fingerprint()
        finger3.add('dummiest', n=5)
        finger3.add('x')
        finger3.add('dummy')

        from fingerprint import FingerprintSet
        fset = FingerprintSet()
        fset += [finger1, finger2, finger3]
        print(f'Jaccard Similarity:{fset.jaccard_sim()}')
        print(f'Bray-Curtis Dissimilarity:{fset.bray_curtis_dis()}')
        print(f'intersect:{fset.jaccard_sim(idx=[0, 3])}')
        exit(0)
