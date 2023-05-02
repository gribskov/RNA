"""=================================================================================================
Convert multi-ct file produced by RNAStructure stochastic to Xios. stochastic samples multiple
structures from a partition function

Michael Gribskov     15 April 2023
================================================================================================="""


# import numpy as np

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
                self.tip_match2(pos)

            pos += 1

        # print out traces for checking
        n = 0
        for t in self.tips:
            n += 1
            print(f'{n}', end='    ')
            while t:
                print(f'{t.pos}:{t.ppos}', end='  ')
                t = t.parent
            print()

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
            for pp in sorted(pairs):
                if pp < pos:
                    continue
                ldif = pos - t[-1][0]
                rdif = t[-1][1] - pp
                # if ldif <= gap and rdif <= gap and rdif > -2:
                if ldif <= gap and abs(rdif) <= gap:
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
                # if a tip is not matched AND pos is more than gap away from the last paired
                # position it can be closed
                remove.append(t)

        # remove un-extensible tips, the stems are still in the list of stems
        for t in remove:
            self.tips.remove(t)

        # check the pairs at this position, any pair not matched to a tip needs a new stem
        last = -gap
        for pp in sorted(pairs, reverse=True):
            found = False
            for m in match:
                if m['pair'] == pp:
                    found = True
            if not found:
                # create new stem
                # last checks to make sure two very near pairs are not both made new tips
                if abs(pp - last) > gap:
                    self.stems.append([(pos, pp)])
                    self.tips.append(self.stems[-1])
                last = pp

        return len(self.tips)

    def tip_match2(self, pos):
        """-----------------------------------------------------------------------------------------
        check if the pos, pairpos can be added to any of the active tips
        if not, create a new stem with the pair,pairpos tuple and make it an active tip

        :param pos: int     position in self.ct
        :return: int        number of active tips
        -----------------------------------------------------------------------------------------"""
        gap = self.gap
        pairs = self.ct[pos]['pair']

        # compare all tips to the pairs at this position and find which pairs can be added to which
        # tips

        for pp in sorted(pairs):
            # examine all the paired positions for pos
            if pp < pos:
                continue
            matched = False
            thisbp = Bp(pos=pos,ppos=pp)
            for t in self.extensible(thisbp):
                if isinstance(t,Bp):
                    # extensible tip
                    self.tips.append(thisbp)
                    thisbp.parent = t
                    if t in self.tips:
                        self.tips.remove(t)
                elif t:
                    break
                else:
                    # if the basepair can't be matched to any tip create a new stem and tip
                    self.stems.append(thisbp)
                    self.tips.append(thisbp)

        # end of loop over basepairs at this position

        return len(self.stems)

    def extensible(self, thisbp):
        """-----------------------------------------------------------------------------------------
        generate positions of all possible extensions using the base-pair pos:ppos

        called once for each paired base so pos is always >= some tip, the latter case occurs when
        a base has two different possible partners that are close to each other

        :param thisbp: bp   new base position to be added to stems
        :return: Bp object  extensible tip
        -----------------------------------------------------------------------------------------"""
        gap = self.gap
        current_tips = self.tips[:]
        matched = False
        while current_tips:
            t = current_tips.pop()
            ldif = thisbp.pos - t.pos
            rdif = t.ppos - thisbp.ppos

            while ldif <= 0 or rdif <= 0:
                # thisbp.ppos must be smaller than t.ppos to be added. If it's bigger, search
                # backwards along the stem to find a position where thisbp.ppos is bigger

                t = t.parent
                if not t:
                    # no more parents - could not match to this stem
                    break

                ldif = thisbp.pos - t.pos
                rdif = t.ppos - thisbp.ppos

            # if t is false there was no extensible basepair on this stem, go on to next tip
            if not t:
                continue

            # otherwise check to see if the possible extension is within gap
            if ldif <= gap and rdif <= gap:
                matched = True
                yield t

        # raise StopIteration
        if not matched:
            yield None
        else:
            yield True

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


def dp(stem):
    """---------------------------------------------------------------------------------------------
    use dp to find the best match of the base-pairs in the stem and generate the vienna string
    :param stem:
    :return:
    ---------------------------------------------------------------------------------------------"""
    import numpy as np

    # find the max and min values in the stem
    left_min, right_max = stem[0]
    left_max, right_min = stem[-1]
    for s in stem:
        left_min = min(left_min, s[0])
        left_max = max(left_max, s[0])
        right_min = min(right_min, s[1])
        right_max = max(right_max, s[1])

    left_len = left_max - left_min + 1
    right_len = right_max - right_min + 1

    # init scoring matrix and fill in paired bases with ones
    score_matrix = np.zeros((left_len, right_len))
    for s in stem:
        score_matrix[s[0] - left_min][right_max - s[1]] = 1

    dir = np.zeros((left_len, right_len))
    for i in range(1, left_len):
        if score_matrix[i][0]:
            dir[i][0] = 0
        else:
            dir[i][0] = 1
    for j in range(1, right_len):
        if score_matrix[0][j]:
            dir[0][j] = 0
        else:
            dir[0][j] = 2

    # Fill in the scoring matrix
    for i in range(1, left_len):
        for j in range(1, right_len):
            match = score_matrix[i - 1][j - 1] + score_matrix[i][j]
            gap1 = score_matrix[i - 1][j]
            gap2 = score_matrix[i][j - 1]
            best = max(match, gap1, gap2)
            if best == match:
                # diagonal
                if score_matrix[i][j]:
                    dir[i][j] = 0
                else:
                    dir[i][j] = 4
            elif best == gap1:
                dir[i][j] = 1
            elif best == gap2:
                dir[i][j] = 2
            score_matrix[i][j] = best

    # Traceback to find the alignment
    vleft = ""
    vright = ""
    i = left_max - left_min
    j = right_max - right_min
    while True:
        # print(f'i:{i}, j:{j}, {score_matrix[i][j]}, {dir[i][j]}')
        if dir[i][j] == 0:
            if score_matrix[i][j]:
                vleft += '('
                vright += ')'
            i -= 1
            j -= 1
        if dir[i][j] == 4:
            vleft += '.'
            vright += '.'
            i -= 1
            j -= 1
        elif dir[i][j] == 1:
            vleft += '.'
            vright += ''
            i -= 1
        elif dir[i][j] == 2:
            vleft += ''
            vright += '.'
            j -= 1

        if i < 0 or j < 0:
            break

        # if score_matrix[i][j] == 1:
        #     vleft += '('
        #     vright += ')'

    # reverse the right vienna string
    vleft = vleft.strip('.')[::-1]
    vright = vright.strip('.')
    # print('\n', score_matrix)
    # print(dir)

    return vleft, vright


class Bp():
    """=============================================================================================
    once base pair in a stem. connected in a linked list where pos < previous_pos and pair_pos >
    previous pair_pos. linked Bp can be branched
    ============================================================================================="""

    def __init__(self, pos=None, ppos=None, ):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.pos = pos
        self.ppos = ppos
        self.parent = None


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
    pos = 0
    for s in struc.ct:
        print(f'{pos:3d}\t{s}')
        pos += 1
    struc.makestems()
    pos = 0
    # for s in struc.stems:
    #     if s[0][0] < s[0][1]:
    #         # print(s)
    #         ave = (s[0][0] + s[-1][0] + s[-1][1] + s[0][1]) / 4.0
    #         print(f'{pos:3d}{ave:6.1f}{s[0][0]:5d}{s[-1][0]:5d}{s[-1][1]:5d}{s[0][1]:5d}', end='  ')
    #         lstem, rstem = dp(s)
    #         print(f'{lstem}     {rstem}')
    #
    #         pos += 1

    exit(0)
