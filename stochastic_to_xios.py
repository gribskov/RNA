"""=================================================================================================
Convert multi-ct file produced by RNAStructure stochastic to Xios. stochastic samples multiple
structures from a partition function

Michael Gribskov     15 April 2023
================================================================================================="""


class Link:
    """=============================================================================================

    Synopsis
        (None yet)

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

    Synopsis
        (None yet)

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
        Read and return one structure from ctfile

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

            if pairpos > 0:
                if pairpos in link.pair:
                    link.pair[pairpos] += 1
                else:
                    link.pair[pairpos] = 1

            if pos == seqlen:
                break

            link = link.next

        return True


# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    ctfilename = 'data/partition.stochastic.ct'
    chain = Chain(ctfilename)

    ct_n = 0
    while True:
        ct_n += 1
        success = chain.ct_read()

        if not success:
            break

    # chain.bp.dump()
    paired = chain.bp.filter(100)
    print(f'{paired} paired positions after filtering')
    chain.bp.dump()

    exit(0)
