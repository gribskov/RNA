"""=================================================================================================
Convert multi-ct file produced by RNAStructure stochastic to Xios. stochastic samples multiple
structures from a partition function

Michael Gribskov     15 April 2023
================================================================================================="""
import copy
import sys
import argparse
import time
from topology import RNAstructure, Stem


class Struc:
    """---------------------------------------------------------------------------------------------
    Struc is used to construct a list of paired positions from the stochastic samples read from the
    ct file. The positions are traversed in order; paired positions are added to stems - the
    stems that can potentially be extended are the tips. At each position all paired bases are
    examined and either added to a a growing stem (tip) or they create a new stem (which is a new
    tip)

    Struc could have been a subclass of topology::Topology, but it seems different enough to
    stand alone
   ----------------------------------------------------------------------------------------------"""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        ctfile      fh, ctfile from RNAStructure Stochastic
        ct          list of positions in sequence, each position has the base and a dictionary of
                    paired positions with counts of the number of ocurrences in the stochastic
                    samples. pair:{pos:count, pos:count, ...}. t[0] is blank because position
                    numbering begins at 1
        tips        list of Bp. used in tracing stems, distinct stems that are currently being
                    extended. Tips are the innermost basepairs in a stem (next to the loop)
        stems       linked list of Bp. contiguous base-paired regions, subject to gap limits
        stemgroups  list of Bp. stem tips that have base pairs in common, each group produces one
                    stem in the final output
        gap         unpaired regions in stems must be <= gap
        id          name of the sequence, from the first line of the ct file
        -----------------------------------------------------------------------------------------"""
        self.ctfile = None
        self.ct = []
        self.tips = []
        self.stems = []
        self.stemgroups = []
        self.gap = 3
        self.id = ''

    def makestems(self):
        """-----------------------------------------------------------------------------------------
        trace stems from the merged ct information. At each position, see if any of the pairs at the
        current position can be added to any of the currently extendable stems (tips)

        :return:
        -----------------------------------------------------------------------------------------"""
        pair = 'pair'
        basepos = 0
        # remember, ct position zero is blank
        for position in self.ct:

            if position[pair]:
                self.tip_match(basepos)

            basepos += 1

        # print out traces for checking
        # for t in self.trace_stem(self.tips[:]):
        #     if t in self.tips:
        #         print(f'start', end='  ')
        #     print(f'{t.pos}:{t.ppos}', end='  ')
        #     if not t.parent:
        #         print(f'end')

        return

    def tip_match(self, pos):
        """-----------------------------------------------------------------------------------------
        check if the pos, pairpos can be added to any of the active tips
        if not, create a new stem with the pair,pairpos tuple and make it an active tip

        :param pos: int     position in self.ct
        :return: int        number of active tips
        -----------------------------------------------------------------------------------------"""
        pairs = self.ct[pos]['pair']

        # compare all tips to the pairs at this position and find which pairs can be added to which
        # tips

        for pp in sorted(pairs):
            # examine all the paired positions for pos
            if pp < pos:
                continue
            # matched = False
            thisbp = Bp(pos=pos, ppos=pp)
            self.stems.append(thisbp)
            for t in self.extensible(thisbp):
                if isinstance(t, Bp):
                    # extensible tip
                    if not thisbp in self.tips:
                        # this basepair not seen before, make it a tip
                        self.tips.append(thisbp)

                    # add parent to basepair (can be multiple) and remove parent from tips
                    thisbp.parent.append(t)
                    if t in self.tips:
                        self.tips.remove(t)
                elif t:
                    # base pair was matched to a tip so you don't have to start a new one
                    break
                else:
                    # if the basepair can't be matched to any tip create a new stem and tip
                    # self.stems.append(thisbp)
                    self.tips.append(thisbp)

        # end of loop over basepairs at this position

        return len(self.stems)

    def final_stems(self, rna, minstem=3):
        """-----------------------------------------------------------------------------------------
        Select the final stem from each possible trace in stemgroups. The selected stem is the first
        stem after sorting by stem length (longer is better) and number of unpaired bases (fewer
        is better)

        :param rna: RNAstructure
        :param minstem:int          minimum number of paired positions in stem
        :return: int                number of stems
        -----------------------------------------------------------------------------------------"""
        stem_n = 0
        for g in self.stemgroups:
            # new = True
            stem_set = []
            stack = []
            for tip in g:
                s = Stem(lbegin=tip.pos, lend=tip.pos, rbegin=tip.ppos, rend=tip.ppos)
                stack.append((s, tip))

            while stack:
                s, pair = stack.pop()
                s.update(pair)

                if pair.parent:
                    # this stem continues
                    for p in pair.parent:
                        stack.append((s.copy(), p))
                else:
                    # this stem is complete
                    stem_set.append(s)

            if not stem_set:
                break
            stems = sorted(stem_set,
                           key=lambda x: (min(x.lend - x.lbegin, x.rend - x.rbegin), -x.unp_n),
                           reverse=True)
            s = stems[0]
            if s.bp_n >= minstem:
                # print(f'{s.formatted()}')
                s.name = stem_n
                stem_n += 1
                rna.stem_list.append(s)

        return stem_n

    def find_groups(self):
        """-----------------------------------------------------------------------------------------
        trace the path from each tip, assigning each basepair to the same group.
        if a trace encounters a basepair in a group, the entire trace from the beginning belongs to
        the older group
        :return:
        -----------------------------------------------------------------------------------------"""
        group = [[]]
        group_n = 0
        for start in self.tips:
            old_group = None
            group[group_n].append(start)
            for t in self.trace_stem([start]):
                if t.group is None:
                    t.group = group_n
                elif t.group != group_n:
                    # this starting point belongs to a previously known group
                    old_group = t.group
                    group[old_group].append(start)
                    group[group_n].remove(start)
                    # group_n -= 1
                    break

            if old_group:
                # need to reset to old groups
                for t in self.trace_stem([start]):
                    t.group = old_group
            else:
                group_n += 1
                group.append([])

        self.stemgroups = group
        return len(self.stemgroups)

    def trace_stem(self, start):
        """-----------------------------------------------------------------------------------------
        generator to trace the linked list from a starting point, considers multiple parents
        :param start:   bp object, beginning position of trace is innermost pair in a stem
        :return:
        -----------------------------------------------------------------------------------------"""
        stack = start
        while stack:
            t = stack.pop()
            yield t

            if t.parent:
                for p in t.parent:
                    stack.append(p)

        return

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

                if t.parent:
                    t = t.parent[0]
                    for tt in range(1, len(t.parent)):
                        current_tips.append(t.parent[tt])
                else:
                    t = None

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
        # if there are comments in the fasta file there will be more tokens, so just look at the first two
        seqlen, self.id = line.rstrip().split()[0:2]
        seqlen = int(seqlen)

        # initialize the list to that tracts the bases and their paired partners
        # elements are a dict of base, pos, [] list of paired positions. the
        # +1 is so that natural coordinates 1 .. seqlen can be used exactly as
        # in the ct file
        if not self.ct:
            self.ct = [{'base': '', 'pair': {}} for _ in range(seqlen + 1)]
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


class Bp:
    """=============================================================================================
    one base pair in a stem. connected in a linked list where pos < previous_pos and pair_pos >
    previous pair_pos. linked Bp can be branched by having multiple parents
    ============================================================================================="""

    def __init__(self, pos=None, ppos=None, ):
        """-----------------------------------------------------------------------------------------
        pos     base on the left side of stem
        ppos    base on the right side of stem
        parent  list of paired positions that can extend this stem
        group   used to gather all intersecting stem traces into a single stem
        -----------------------------------------------------------------------------------------"""
        self.pos = pos
        self.ppos = ppos
        self.parent = []
        self.group = None


class Stem(Stem):
    """=============================================================================================
    Subclasses topology.py:Stem
    holds the entire description of a stem
    ============================================================================================="""

    def __init__(self, lbegin=0, lend=0, rbegin=0, rend=0):
        """-----------------------------------------------------------------------------------------
        lbegin, lend        begin and end position of left halfstem
        rbegin, rend        begin and end position of right halfstem
        bp_n                number of paired position
        unp_n               number of unpaired positions, counted on both halfstems
        -----------------------------------------------------------------------------------------"""
        super(Stem, self).__init__()
        self.lbegin = lbegin
        self.lend = lend
        self.rbegin = rbegin
        self.rend = rend
        self.bp_n = 0
        self.unp_n = 0

    def copy(self):
        return copy.copy(self)

    def update(self, pair):
        """-----------------------------------------------------------------------------------------
        add the current pair to the stem, adding spaces as necessary to the vienna strings
        :param pair:
        :return:
        -----------------------------------------------------------------------------------------"""
        self.bp_n += 1
        for p in range(pair.pos, self.lbegin - 1):
            self.lvienna += '.'
            self.unp_n += 1
        self.lvienna += '('
        self.lbegin = pair.pos

        for p in range(self.rend, pair.ppos - 1):
            self.rvienna += '.'
            self.unp_n += 1
        self.rvienna += ')'
        self.rend = pair.ppos

        return


# --------------------------------------------------------------------------------------------------
# main program functions
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    def arg_formatter(prog):
        """-----------------------------------------------------------------------------------------
        Set up formatting for help. Used in arg_get.

        :param prog:
        :return: argparse formatter class
        -----------------------------------------------------------------------------------------"""
        return argparse.HelpFormatter(prog, max_help_position=60, width=120)


    def getoptions():
        """-----------------------------------------------------------------------------------------
        Command line arguments:
        project: project directory, all stages will be in subdirectories

        r, restart      remove previous results and perform all stages
        c, continue     continue existing run, skipping commands that are marked as complete

        options:
        j, jobs         number of concurrent jobs to run
        l, log          directory for log files

        :return: command line namespace
        -----------------------------------------------------------------------------------------"""
        commandline = argparse.ArgumentParser(
            description='Calculate XIOS from RNAStructure::Stochastic',
            formatter_class=arg_formatter)

        commandline.add_argument('input_ct', type=str,
                                 help='CT file produced by Stochastic.exe (%(default)s)',
                                 default='rna.ct')

        commandline.add_argument('output_xios', type=str,
                                 help='XIOS file (%(default)s)',
                                 default='rna.xios')

        commandline.add_argument('-m', '--minstem',
                                 help='Minimum number of paired bases in a stem (%(default)s)',
                                 default=3, type=int)

        args = commandline.parse_args()
        # return vars(args)  # convert namespace to dict
        return args


    # ----------------------------------------------------------------------------------------------
    # actual main program starts here
    # ----------------------------------------------------------------------------------------------
    opt = getoptions()
    now = time.localtime()
    sys.stdout.write(f'stochastic_to_xios {time.asctime(now)}\n')
    sys.stdout.write(f'\tinput ct file: {opt.input_ct}\n')
    sys.stdout.write(f'\toutput xios file: {opt.output_xios}\n')
    sys.stdout.write(f'\tminimum stem size: {opt.minstem}\n')

    struc = Struc()
    struc.ctfile = open(opt.input_ct, 'r')
    struc.ct_read_all()
    struc.filter(56)

    rna = RNAstructure()
    rna.sequence_id = struc.id
    comment = f'creation_date: {time.asctime(now)}\n'
    comment += f'Program: xios_from_stochastic\n'
    comment += f'input_file: {opt.input_ct}\n'
    comment += f'minimum stem: {opt.minstem}'
    rna.comment.append(comment)

    struc.makestems()
    group_n = struc.find_groups()
    sys.stdout.write(f'\n\t{group_n} groups found\n')
    stem_n = struc.final_stems(rna, opt.minstem)
    sys.stdout.write(f'\t{stem_n} stems selected\n\n')
    rna.adjacency_from_stemlist()
    rna.edgelist_from_adjacency(include="ijo", whole=False)
    rna.XIOSwrite(open(opt.output_xios, 'w'))

    exit(0)
