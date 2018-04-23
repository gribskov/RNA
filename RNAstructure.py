class RNAstructure:
    """=============================================================================================
    RNA structure
    a single RNA structure
    what is a structure?  well might you ask.
    ============================================================================================="""

    def __init__(self):
        self.sequence = ''
        self.pair = []  # base number of the paired base
        self.stemlist = []
        self.length = 0
        self.id = None

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
        with open(filename, 'r') as ct:
            line = ct.readline()
            # print('firstline:', line)
            # print('field:',field)

            if line.find('ENERGY') >= 0:
                # TODO: not sure what this is checking
                pass
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

        if instem:
            # if you  end in a stem, save it before closing
            stem.trimVienna()
            self.stemlist.append(stem)

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

    exit(0)

""" correct by manual inspection
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
