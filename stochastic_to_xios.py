"""=================================================================================================
Convert multi-ct file produced by RNAStructure stochastic to Xios. stochastic samples multiple
structures from a partition function

Michael Gribskov     15 April 2023
================================================================================================="""


class Struc:
    """---------------------------------------------------------------------------------------------
    Struc is used to construct a list of paired postions from the stochastic samples read from the
    ct file. The positions are traversed in order; paired positions are added to stems - the
    stems that can potentially be extended are the tips. At each position all paired bases are
    examined and either added to a a growing stem (tip) or they create a new stem (which is a new
    tip)
   ----------------------------------------------------------------------------------------------"""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        ctfile  fh, ctfile from RNAStructure Stochastic
        ct      list of positions in sequence, each position has the base and a dictionary of paired
                positions and a count of number of ocurrences pair:{pos:count, pos:count, ...}
                ct[0] is blank because position numbering begins at 1
        tips    used in tracing stems, distinct stems that are currently being extended
        stems   contiguous base-paired regions, subject to gap
        gap     stems must have fewer than gap contiguous unpaired positions ((((..((....))..))..))
                is ok
        -----------------------------------------------------------------------------------------"""
        self.ctfile = None
        self.ct = []
        self.tips = []
        self.stems = []
        self.gap = 3

    def makestems(self):
        """-----------------------------------------------------------------------------------------
        trace stems from the merged ct information. At each position, see if any of the pairs at the
        current position can be added to any of the currently extendable stems (tips)

        :return:
        -----------------------------------------------------------------------------------------"""
        pair = 'pair'
        pos = 0
        # remember, ct position zero is blank
        for position in self.ct:

            if position[pair]:
                self.tip_match(pos)

            pos += 1

        return

    def tip_update(self, pos, pairpos, gap=3):
        """-----------------------------------------------------------------------------------------
        Tips are the contguous base-paired regions (stems) that are currently available by
        extension. Tips cease to be extensible when a unpaired region (gap) of gap or more is
        encountered.

        remove any active tips where adding pair, pairpos would create a bubble of gap or more on
        both sides of the stem

        :param gap: int     contiguous unpaired regions must < gap long
        :return: int        current number of tips
        -----------------------------------------------------------------------------------------"""
        if not self.tips:
            # no active tips create a new stem and return
            self.stem_new(pos, pairpos)
            return len(self.tips)

        match = False
        self.match_find(pos)

        return len(self.tips)

    def stem_new(self, pos, pairpos):
        """-----------------------------------------------------------------------------------------
        Stem is a list of positions and paired position that can be joined subject to the limit on
        gaps
        :param pos: int
        :param pairpos: int
        :return: list           last paired position in stem
        -----------------------------------------------------------------------------------------"""
        self.stems.append([(pos, pairpos)])
        self.tips.append(self.stems[-1])
        return self.stems[-1]

    def tip_match(self, pos):
        """-----------------------------------------------------------------------------------------
        check if the pos, pairpos can be added to any of the active tips
        if not, create a new stem with the pair,pairpos tuple and make it an active tip

        :param pos: int     position in self.ct
        :return: int        number of active tips
        -----------------------------------------------------------------------------------------"""
        gap = self.gap
        pairs = self.ct[pos]['pair']
        match = []

        # compare all tips to the pairs at this position and find which pairs can be added to which
        # tips
        for t in self.tips:
            for pp in pairs:
                ldif = pos - t[-1][0]
                rdif = t[-1][1] - pp

                # if pos - t[-1][0] <= gap and pp - t[-1][1] <= gap and t[-1][1] - pp <= gap:
                if ldif <= gap and rdif <= gap and ldif > 0 and rdif > 0:
                    match.append({'pos': pos, 'pair': pp, 'tip': t})

        # make a list of tips that will be removed because they are no longer extensible
        remove = []
        for t in self.tips:
            found = False
            for m in match:
                if m['tip'] == t:
                    t.append((m['pos'], m['pair']))
                    found = True

            if not found and pos - t[-1][0] > gap:
                # if a tip is not matched AND is pos is more than gap away from the last paired
                # position it can be closed
                remove.append(t)

        # remove un-extensible tips, the stems are still in the list of stems
        for t in remove:
            self.tips.remove(t)

        # check the pairs at this position, any pair not matched to a tip needs a new stem
        for pp in pairs:
            found = False
            for m in match:
                if m['pair'] == pp:
                    found = True
            if not found:
                # create new stem
                self.stems.append([(pos, pp)])
                self.tips.append(self.stems[-1])

        return len(self.tips)

    def ct_read_all(self):
        """-----------------------------------------------------------------------------------------
        Read and return all structures from ctfile

        :return:
        -----------------------------------------------------------------------------------------"""
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
        """-----------------------------------------------------------------------------------------
        remove pairs whose count < n

        :param n: int   minimum count for base pairs
        :return: int    number of positions with pairs
        -----------------------------------------------------------------------------------------"""
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


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
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
            print(s)
            ave = (s[0][0] + s[-1][0] + s[-1][1] + s[0][1]) / 4.0
            print(f'{pos:3d}{ave:6.1f}{s[0][0]:5d}{s[-1][0]:5d}{s[-1][1]:5d}{s[0][1]:5d}', end='  ')

            # construct vienna string
            left = right = ''
            lpos = s[0][0]
            rpos = s[0][1]
            for pos, pairpos in s:
                while lpos < pos:
                    left += '.'
                    lpos += 1
                while rpos > pairpos:
                    right += '.'
                    rpos -= 1
                if pos == lpos:
                    left += '('
                    lpos += 1
                if pairpos == rpos:
                    right += ')'
                    rpos -= 1

            print(f'{left}    {right[::-1]}')

            pos += 1

    exit(0)
