"""=================================================================================================
Convert multi-ct file produced by RNAStructure stochastic to Xios. stochastic samples multiple
structures from a partition function

Michael Gribskov     15 April 2023
================================================================================================="""


class Link:
    """=============================================================================================
    A link are is two paired bases, basically a linked list of pairs

    Michael Gribskov
    ============================================================================================="""

    def __init__(self, n=None):
        """-----------------------------------------------------------------------------------------
         link constructor
        -----------------------------------------------------------------------------------------"""
        self.base = ''
        self.pos = n
        self.pair = {}
        self.next = None

    @classmethod
    def build_chain(cls, len):
        """-----------------------------------------------------------------------------------------
        Build a linked list of the desired length. because RNA sequences are one based, the first
        link is position 1 and the last is position len

        :param len: int     length of chain = length of RNA sequence
        :return:
        -----------------------------------------------------------------------------------------"""
        root = Link()
        l = root
        for pos in range(1, len + 1):
            l.next = Link(n=pos)
            l = l.next

        return root.next

    def dump(self):
        """-----------------------------------------------------------------------------------------
        print out the linked list

        :return: int    number of links
        -----------------------------------------------------------------------------------------"""
        link_n = 0
        link = self
        while True:
            link_n += 1
            print(f'{link.pos:3d}\t{link.base}')
            for pp in link.pair:
                print(f'\t{pp}\t{link.pair[pp]}')

            if not link.next:
                break
            link = link.next

        return link_n

    def filter(self, n):
        """-----------------------------------------------------------------------------------------
        remove pairs whose count < n

        :param n: int   minimum count for base pairs
        :return: int    number of positions with pairs
        -----------------------------------------------------------------------------------------"""
        link_n = 0
        link = self
        pair_n = 0
        while True:
            remove = []
            for pp in link.pair:
                if link.pair[pp] < n:
                    remove.append(pp)

            for pp in remove:
                del (link.pair[pp])

            if link.pair:
                pair_n += 1

            if not link.next:
                break

            link = link.next

        return pair_n


class Chain:
    """=============================================================================================


    Michael Gribskov
    ============================================================================================="""

    def __init__(self, ctfilename=''):
        """-----------------------------------------------------------------------------------------
        chain constructor
        -----------------------------------------------------------------------------------------"""
        self.id = ''
        self.ctfilename = ctfilename
        self.ctfile = None
        self.bp = None

        if self.ctfilename:
            self.ctfile = open(ctfilename, 'r')

    def ct_read(self):
        """-----------------------------------------------------------------------------------------
        Read and return the next ct structure in the file, example of format below.

          360  rnasep_a2.Neisseria_meningitidis.fa
        1 C       0    2    0    1
        2 G       1    3    0    2
        3 G       2    4  354    3
        4 G       3    5  353    4
        5 A       4    6  352    5
        6 C       5    7  351    6
        from https://rna.urmc.rochester.edu/Text/File_Formats.html#CT
        line 0:
            Start of first line: number of bases in the sequence
            End of first line: title of the structure
        Each of the following lines, one/base
            Base number: index n
            Base (A, C, G, T, U, X)
            Index n-1
            Index n+1
            Number of the base to which n is paired. No pairing is indicated by 0 (zero).
            Natural numbering. RNAstructure ignores the actual value given in natural numbering, so it is easiest to
            repeat n here.

        :param ct: file     ct file open for reading
        :return: dict       ct information

        :return:
        -----------------------------------------------------------------------------------------"""
        ctfile = self.ctfile

        # first line
        line = ctfile.readline()
        try:
            seqlen, self.id = line.rstrip().split()
        except ValueError:
            # end of file
            return False

        seqlen = int(seqlen)

        if not self.bp:
            # initialize paired sequence linked list
            self.bp = Link.build_chain(seqlen)

        link = self.bp
        while True:
            line = ctfile.readline()
            pos, base, x0, x1, pairpos, x2 = line.rstrip().split()
            pos = int(pos)
            pairpos = int(pairpos)
            link.base = base

            # to get only stems where left < right use
            # if pairpos > pos:
            if pairpos > 0:
                if pairpos in link.pair:
                    link.pair[pairpos] += 1
                else:
                    link.pair[pairpos] = 1

            if pos == seqlen:
                break

            link = link.next

        return True


def trace(link, gap=3):
    """-----------------------------------------------------------------------------------------------------------------
    try to trace the stems
    :param link:
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    lbeg = 'lbeg'
    lend = 'lend'
    rbeg = 'rbeg'
    rend = 'rend'
    stack = [(link, 0)]
    stem = {}
    while stack:
        link, pp = stack.pop()

        # skip unpaired
        if not link.pair:
            if link.next:
                stack.append((link.next, 0))
            print(f'{link.pos:3d}\t{link.base}\tunpaired')
        else:
            if stem:
                # there is an existing stem
                if stem[lend] + gap >= link.pos and stem[rbeg] - gap <= link.pair[pp]:
                    stem[lend] = link.pos
                    stem[rbeg] = link.pair[pp]
                else:
                    # start new stem
                    print(stem)
                    stem[lbeg] = stem[lend] = link.pos
                    stem[rbeg] = stem[rend] = link.pair[pp]
            if link.next:
                stack.append((link.next, 0))

    return

    # TODO probably better to actually condense the list


class Struc:
    """-----------------------------------------------------------------------------------------------------------------

   ----------------------------------------------------------------------------------------------------------------- """

    def __init__(self):
        self.ctfile = None
        self.ct = []
        self.tips = []
        self.stems = []
        self.gap = 3

    def ct_read(self):
        """-----------------------------------------------------------------------------------------
        Read and return one structure from ctfile

        :return:
        -----------------------------------------------------------------------------------------"""
        ctfile = self.ctfile

        # first line
        line = ctfile.readline()
        try:
            seqlen, self.id = line.rstrip().split()
            seqlen = int(seqlen)
        except ValueError:
            # end of file
            return False

        while True:
            line = ctfile.readline()
            pos, base, x0, x1, pairpos, x2 = line.rstrip().split()
            pos = int(pos)
            pairpos = int(pairpos)

            # to get only stems where left < right use
            # if pairpos > pos:
            if pairpos > pos:
                active_n = self.tip_update(pos, pairpos)

            if pos == seqlen:
                break

        return True

    def makestems(self):
        """-------------------------------------------------------------------------------------------------------------
        trace stems from the merged ct files

        :return:
        -------------------------------------------------------------------------------------------------------------"""
        pair = 'pair'
        pos = 0
        # remember, ct position zero is blank
        for position in self.ct:

            if position[pair]:
                self.tip_match(pos)

            pos += 1
            # for pairpos in position[pair]:
            #     active_n = self.tip_update(pos, pairpos)

    def tip_update(self, pos, pairpos, gap=3):
        """-------------------------------------------------------------------------------------------------------------
        remove any active tips where adding pair, pairpos would create a bubble of gap or more on
        both sides of the stem

        :param gap:
        :return:
        -------------------------------------------------------------------------------------------------------------"""
        if not self.tips:
            # no active tips create a new stem and return
            self.stem_new(pos, pairpos)
            return len(self.tips)

        match = False
        self.match_find(pos)
        # for s in self.tips:
        #     # compare last position to pair, pairpos
        #     send = s[-1]
        #     if pos - send[0] > gap or pairpos - send[1] > gap or send[1] - pairpos > gap:
        #         # gap is too big
        #         self.tips.remove(s)
        #     else:
        #         match = True
        #         s.append((pos, pairpos))
        #
        # if not match:
        #     # no match to any active tip, create a new stem
        #     self.stem_new(pos, pairpos)

        return len(self.tips)

    def stem_new(self, pos, pairpos):
        """-------------------------------------------------------------------------------------------------------------

        :param pos:
        :param pairpos:
        :return:
        -------------------------------------------------------------------------------------------------------------"""
        self.stems.append([(pos, pairpos)])
        self.tips.append(self.stems[-1])
        return self.stems[-1]

    def tip_match(self, pos):
        """-------------------------------------------------------------------------------------------------------------
        check if the pos, pairpos can be added to any of the active tips
        if not, create a new stem with the pair,pairpos tuple and make it an
        active tip

        :param pos:
        :param pairpos:
        :return:
        -------------------------------------------------------------------------------------------------------------"""
        gap = self.gap
        pairs = self.ct[pos]['pair']
        # any_match = False
        match = []
        for t in self.tips:
            for pp in pairs:
                if pos - t[-1][0] <= gap and pp - t[-1][1] <= gap and t[-1][1] - pp <= gap:
                    match.append({'pos': pos, 'pair': pp, 'tip': t})
                    # any_match = True

        remove = []
        for t in self.tips:
            found = False
            for m in match:
                if m['tip'] == t:
                    # TODO this has to check the position and not terminate until gap bases away
                    t.append((m['pos'], m['pair']))
                    found = True

            if not found and pos - t[-1][0] > gap:
                # if a tip is not matched AND is pos is more than gap it can be closed
                remove.append(t)

        for t in remove:
            self.tips.remove(t)

        # check the pairs, any pair not found needs a new stem
        for pp in pairs:
            found = False
            for m in match:
                if m['pair'] == pp:
                    found = True
            if not found:
                # create new stem
                self.stems.append([(pos, pp)])
                self.tips.append(self.stems[-1])

        return

    def ct_read_all(self):
        """-------------------------------------------------------------------------------------------------------------
        Read and return all structures from ctfile

        :return:
        -------------------------------------------------------------------------------------------------------------"""
        ctfile = self.ctfile

        # first line
        line = ctfile.readline()
        seqlen, self.id = line.rstrip().split()
        seqlen = int(seqlen)

        # initialize the list to that tracts the bases and their paired partners
        # elements are a dict of base, pos, [] list of paired positions. the
        # +1 is so that natural coordinates 1 .. seqlen can be used exactly as
        # in the ct file
        if not self.ct:
            self.ct = [{'base': '', 'pair': {}} for i in range(seqlen + 1)]
        base = 'base'
        pair = 'pair'
        while True:
            line = ctfile.readline()
            if not line:
                # eof
                break

            # skip lines with length and title (only two tokens, 6 expected)
            try:
                pos, nuc, x0, x1, pairpos, x2 = line.rstrip().split()
            except ValueError:
                continue

            pos = int(pos)
            pairpos = int(pairpos)
            thispos = self.ct[pos]
            thispos[base] = nuc

            # to get only stems where left < right use: if pairpos > pos:
            if pairpos > 0:
                if pairpos in thispos[pair]:
                    thispos[pair][pairpos] += 1
                else:
                    thispos[pair][pairpos] = 1

        return True

    def filter(self, n):
        """-------------------------------------------------------------------------------------------------------------
        remove pairs whose count < n

        :param n: int   minimum count for base pairs
        :return: int    number of positions with pairs
        -------------------------------------------------------------------------------------------------------------"""
        pair_n = 0
        for pos in range(len(self.ct)):
            remove = []
            pairs = self.ct[pos]['pair']
            if pairs:
                pair_n += 1
            for pp in pairs:
                if pairs[pp] < n:
                    remove.append(pp)

            for pp in remove:
                del pairs[pp]

        return pair_n


# ----------------------------------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    struc = Struc()
    ctfilename = 'data/partition.stochastic.ct'
    # ctfilename = 'data/stochastic.330.ct'
    struc.ctfile = open(ctfilename, 'r')
    struc.ct_read_all()
    struc.filter(56)
    pos = 1
    for s in struc.ct:
        print(f'{pos:3d}\t{s}')
        pos += 1
    struc.makestems()
    pos = 0
    for s in struc.stems:
        if len(s) > 3 and s[0][0] < s[0][1]:
            ave = (s[0][0] + s[-1][0] + s[-1][1] + s[0][1]) / 4.0
            print(f'{pos:3d}{ave:6.1f}{s[0][0]:5d}{s[-1][0]:5d}{s[-1][1]:5d}{s[0][1]:5d}')
            pos += 1
    # chain = Chain(ctfilename)
    #
    # ct_n = 0
    # while True:
    #     ct_n += 1
    #     success = chain.ct_read()
    #
    #     if not success:
    #         break
    #
    # # chain.bp.dump()
    # paired = chain.bp.filter(100)
    # print(f'{paired} paired positions after filtering')
    # chain.bp.dump()
    # trace(chain.bp)

    exit(0)
