class RNAstructure:
    """=============================================================================================
    RNA structure
    a single RNA structure
    what is a structure?  well might you ask.
    ============================================================================================="""

    def __init__(self):
        self.sequence = ''
        self.length = 0
        self.energy = None  # not defined for probknot
        self.pair = []  # base number of the paired base
        self.stemlist = []
        self.adjacency = None
        self.id = None
        self.filename = None

    def CTRead(self, filename):
        """-----------------------------------------------------------------------------------------
        Read in RNAstructure CT file

        format example:
          240  ENERGY = -49.3   mr_s129
            1 A       0    2    0    1
            2 A       1    3    0    2
            3 C       2    4    0    3
            4 C       3    5    0    4
            5 A       4    6    0    5

            base_number base previous_base next_base paired_base base_number

        usage
            rna.CTRead(filename)
        :param filename: string, filename with CT formatted structure
        :return: integer, number of bases
        -----------------------------------------------------------------------------------------"""
        nbase = 0
        self.filename = filename  # TODO should strip directory path from filename
        with open(filename, 'r') as ct:
            line = ct.readline()
            # print('firstline:', line)
            # print('field:',field)

            if line.find('ENERGY') >= 0:
                field = line.split()
                self.energy = field[3]

            else:
                # probknot file
                field = line.split()
                self.length = int(field[0])
                self.id = field[1]

            self.pair = [0] * (self.length + 1)

            for line in ct:
                n, base, prev, following, pair, n2 = line.split()
                # print('n:', n, 'base:', base, 'pref:', prev, 'next:', next, 'pair:', pair, 'n2:',
                #       n2)
                self.sequence += base
                if pair != '0':
                    self.pair[int(pair)] = int(n)
                    self.pair[int(n)] = int(pair)
                nbase += 1

        return nbase

    def __str__(self):
        """-----------------------------------------------------------------------------------------
        string representation of a structure
        :return: string
        -----------------------------------------------------------------------------------------"""
        rnastr = 'RNA:\n'
        for key in self.__dict__:
            rnastr += '{0} = {1}\n'.format(key, self.__dict__[key])

        return rnastr

    def stemListGet(self, unpaired=2):
        """-----------------------------------------------------------------------------------------
        Construct the stemlist from the paired base list in stem.pair
        :return: integer, number of stems in stemlist
        -----------------------------------------------------------------------------------------"""
        maxgap = unpaired + 1
        nstem = 0
        instem = False
        for pos in range(0, len(self.pair) - 1):
            if self.pair[pos] == 0 or self.pair[pos] < pos:
                continue

            if instem:
                # currently in a stem
                lgap = pos - stem.lend - 1
                rgap = stem.rbegin - self.pair[pos] - 1
                if lgap >= maxgap or rgap >= maxgap:
                    # gap is too big, end old stem
                    stem.trimVienna()
                    instem = False
                else:
                    # extend current stem
                    stem.lend = pos
                    stem.rbegin = self.pair[pos]
                    stem.lvienna += '{}('.format('.' * lgap)
                    stem.rvienna = '){}'.format('.' * rgap) + stem.rvienna
                    continue

            # not in a stem, start a new stem
            stem = Stem()
            nstem += 1
            self.stemlist.append(stem)
            instem = True

            stem.lbegin = pos
            stem.rend = self.pair[pos]
            stem.lend = pos
            stem.rbegin = self.pair[pos]
            stem.lvienna = '('
            stem.rvienna = ')'

        # if instem:
        # if you  end in a stem, clean up the Vienna string
        # stem.trimVienna()
        # self.stemlist.append(stem)

        return nstem

    def stemlistFormat(self):
        """-----------------------------------------------------------------------------------------
        Returns a string with the stemlist formatted according to Stem.formatted()
        :return: string
        -----------------------------------------------------------------------------------------"""
        n = 0
        stemstr = ''
        for stem in self.stemlist:
            n += 1
            stemstr += '{0}\t{1}\n'.format(n, stem.formatted())
        return stemstr

    def adjacencyGet(self):
        """-----------------------------------------------------------------------------------------
        Calculate an adjacency matrix from a stemlist. assumes stesm are ordered by the beginning
        of the left half-stem.
        :return: dict, keys are edge types, values are counts
        -----------------------------------------------------------------------------------------"""
        nstem = len(self.stemlist)

        edges = {'i': 0, 'j': 0, 'o': 0, 's': 0, 'x': 0}
        a = [[0 for _ in range(nstem)] for _ in range(nstem)]

        for i in range(nstem):
            stem_i = self.stemlist[i]

            for j in range(i + 1, nstem):
                stem_j = self.stemlist[j]

                if stem_i.rend < stem_j.lbegin:
                    # serial edge
                    a[i][j] = 's'
                    a[j][i] = 's'
                    edges['s'] += 1

                elif stem_i.lend < stem_j.lbegin and stem_i.rend < stem_j.rbegin:
                    # overlap edge (pseudoknot)
                    a[i][j] = 'o'
                    a[j][i] = 'o'
                    edges['o'] += 1

                elif stem_i.lend < stem_j.lbegin and stem_i.rbegin > stem_j.rend:
                    # included edge (directed, j is inside i)
                    a[i][j] = 'i'
                    a[j][i] = 'j'
                    edges['i'] += 1

                else:
                    # excluded edge (stems overlap)
                    a[i][j] = 'x'
                    a[j][i] = 'x'
                    edges['x'] += 1

        self.adjacency = a
        return edges

    def adjacencyFormat(self):
        """-----------------------------------------------------------------------------------------

        :return: string, formatted version of adjacency matrix
        -----------------------------------------------------------------------------------------"""
        adjstr = ''
        width = 4  # TODO set based on matrix size?
        if not self.adjacency:
            adjstr = 'None'
            return adjstr

        size = len(self.adjacency)
        for i in range(size):
            for j in range(size):
                # must cast to str so zeros format correctly
                adjstr += '{:{}}'.format(str(self.adjacency[i][j]), width)
            adjstr += '\n'

        return adjstr

    def edgelist(self, include="ijo", whole=False):
        """-----------------------------------------------------------------------------------------
        An edgelist is an array of lists.  each row corresponds to a stem (vertex).  The values are
        tuples with the number and type of nodes with edges
        :return: list, edgelist
        -----------------------------------------------------------------------------------------"""
        elist = []
        if not self.adjacency:
            return elist

        size = len(self.adjacency)
        a = self.adjacency
        for i in range(size):
            e = []
            elist.append(e)
            begin = i + 1
            if whole:
                begin = 0

            for j in range(begin, size):
                if i == j:
                    continue
                if a[i][j] in include:
                    e.append([j, a[i][j]])
        return elist


class Stem:
    """=============================================================================================
    coordinates of stem regions (including allowed gaps) and the corresponding Vienna strings.
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Stem constructor
        begin-end coordinates are small->big for left side
        end-begin for right side
        -----------------------------------------------------------------------------------------"""
        self.lbegin = 0
        self.lend = 0
        self.rbegin = 0
        self.rend = 0
        self.lvienna = ''
        self.rvienna = ''

    def formatted(self):
        """-----------------------------------------------------------------------------------------
        Return a formatted version of the stem:
            left begin pos
            left end pos
            right begin pso
            right end pos
            left Vienna string
            right Vienna string
        :return: string, with stem coordinates and Vienna structure
        -----------------------------------------------------------------------------------------"""
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(self.lbegin, self.lend, self.rbegin, self.rend,
                                                     self.lvienna, self.rvienna)

    def trimVienna(self):
        """-----------------------------------------------------------------------------------------
        Removes leading and trailing unpaired positions (.) from vienna strings
        :return: True
        -----------------------------------------------------------------------------------------"""
        self.lvienna = self.lvienna.rstrip('.')
        self.rvienna = self.rvienna.lstrip('.')
        return True


# ==================================================================================================
# main/test
# ==================================================================================================

if __name__ == '__main__':
    rna = RNAstructure()
    rna.CTRead('data/mr_s129.probknot.ct')
    rna.stemListGet()
    # print(rna)
    print('Stemlist\n')
    print(rna.stemlistFormat())

    edges = rna.adjacencyGet()
    print('edges', edges)
    print('\nAdjacency matrix\n')
    print(rna.adjacencyFormat())

    print('\nEdgelist\n')
    e = rna.edgelist()
    for i in range(len(e)):
        print('{}:\t{}'.format(i, e[i]))

    exit(0)

""" correct structure of mr_s129 by manual inspection
29-31:231-229
33-37:225-221
41-44:220-217
48-52:211-207
55-58:68-65
65-68:58-55
76-80:192-188
84-91:184-177
95-98:155-152
100-103:150-147
106-108:145-143
117-119:130-128
129-130:119-117
143-145:108-106
147-150:103-100
152-155:98-95
163-165:176-174
174-176:165-163
177-184:91-84
188-192:80-76
201-205:216-212
207-211:52-48
212-216:205-201
217-225:44-33
229-231:31-29
"""
