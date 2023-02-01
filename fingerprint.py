import sys
import json
import datetime
import yaml


# note install PyYAML


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
        mm = self.motif_matrix
        j_sim = [[0 for _ in range(len(self))] for _ in range(len(self))]

        for i in range(0,len(self)):
            for j in range(i+1, len(self)):
                intersect = 0
                union = 0
                vi = mm[i]
                vj = mm[j]
                for k in range(0,len(self.i2motif)):
                    intersect += vi[k] and vj[k]
                    union += vi[k] or vj[k]

                try:
                    j_sim[i][j] = intersect / union
                except ZeroDivisionError:
                    j = 0

        return


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

    def select(self, selected):
        """-----------------------------------------------------------------------------------------
        Filter each fingeprint in the fingerprint set retaining only those in selected. If selected
        is empty, do nothing and return

        :param selected: list       motif names
        :return:
        -----------------------------------------------------------------------------------------"""
        if not selected:
            return

        for fpt in self:
            ok_motifs = {}
            for motif in fpt.motif:
                if motif in selected:
                    ok_motifs[motif] = fpt.motif[motif]
                else:
                    fpt.count -= 1

            fpt.motif = ok_motifs

        return

    def select_binary(self, selected=None):
        """-----------------------------------------------------------------------------------------
        Filter each fingeprint in the fingerprint set retaining only those in selected. If selected
        is empty, do nothing and return

        :param selected: list       motif names
        :return: True
        -----------------------------------------------------------------------------------------"""
        # set up translation between motif id and numeric index
        if not selected:
            # all motifs
            motif_n = self.index_motifs()
        else:
            # selected motifs
            self.i2motif = selected
            for i in range(len(selected)):
                self.motif2i[selected[i]] = i
            motif_n = len(selected)

        motif2i = self.motif2i
        i2motif = self.i2motif
        motif_matrix = []
        for fpt in self:
            motif_vector = [0 for _ in range(motif_n)]
            motif_matrix.append(motif_vector)
            for motif in fpt.motif:
                if motif in motif2i:
                    motif_vector[motif2i[motif]] = 1
            print(f'sum {sum(motif_vector)}')

        self.motif_matrix = motif_matrix
        return True

    def index_motifs(self):
        """-----------------------------------------------------------------------------------------
        make an index of all the motifs found in the fingerprintset.
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

    fset = FingerprintSet()
    fset += [finger1, finger2, finger3]
    print(f'Jaccard Similaity:{fset.jaccard_sim()}')
    print(f'Bray-Curtis Dissimilarity:{fset.bray_curtis_dis()}')
    print(f'intersect:{fset.jaccard_sim(idx=[0, 3])}')
    exit(0)
