from topology import Topology, PairRNA
from xios import Gspan, MotifDB


def check_uniqueness(unique, topo):
    """---------------------------------------------------------------------------------------------
    some duplicates can be distinguished simply by checking the string representation of the pairs.
    check in forward and reverse direction since XIOS does not distinguish orientation. unique is
    simple a dictionary with keys corresponding to the already seen pair topologies stored in both
    forward and reverse directions

    :param unique: dict of PairRNA
    :param topo: PairRNA
    :return: boolean, True if structure appears unique (check further with Gspan)
    ---------------------------------------------------------------------------------------------"""
    s = str(topo)
    r = str(topo.reverse())
    if s in unique:
        # print(f'\tnotunique {s}')
        unique[s] += 1
        unique[r] += 1
        return False

    # if unseen store both forward and reverse, this minimizes the number of reversals
    unique[s] = 1
    unique[r] = 1
    topo.reverse()  # change back to original orientation

    return True


def all_pairs(n):
    """---------------------------------------------------------------------------------------------
    Generator that yields all connected graphs with n stems as PairRNA objects.

    :param n: int, number of stems
    :return:
    ---------------------------------------------------------------------------------------------"""
    maxpos = n * 2

    # available locations in descending order
    available = [maxpos - x - 1 for x in range(maxpos)]

    # each item on the stack is a list of pairs and a list of the still available positions
    # initialize stack with empty list and all locations
    stack = [[[], available]]

    while stack:
        pairs, avail = stack.pop()
        # the left location in the next pair is always the next lowest location
        # for the first pair, left is always 0
        left = avail.pop()

        rp = len(avail) - 1
        # rp = 0
        # for rp in range(len(avail)):
        while rp >= 0:
            right = avail[rp]
            rp -= 1
            # all other locations are available for the right position, push all on the stack
            topo = pairs[:] + [left, right]
            newpairs = PairRNA(topology=topo)
            # newavail = avail[:].remove(right)
            newavail = avail[:]
            newavail.remove(right)
            if newavail:
                # check if the graph is/will be disconnected, the graph is disconnected all
                # positions less than the current length have been used,
                # i.e. avail[-1] == next available position
                if len(newpairs.pairs) == newavail[-1]:
                    # print(f'\tdisconnected {newpairs})')
                    continue

                else:
                    # connected, push on stack
                    stack.append([newpairs.pairs, newavail])
            else:
                # reached the desired size, don't push
                yield newpairs

    return


# ===================================================================================================
# main
# ===================================================================================================
if __name__ == '__main__':

    pair_idx = {}  # index of unique pair graphs

    level_max = 7
    motif = MotifDB()

    count = [0, 0]
    cumul = 0

    for i in range(2, level_max + 1):
        count.append(0)
        possible = 0
        assymetric = 0

        for p in all_pairs(i):
            possible += 1
            rp = p.reverse()
            if rp.pairs > p.pairs:
                # > test makes sure that symmetric structures get processed
                # print(i, p, rp, 'R')
                continue
            else:
                assymetric += 1
                # print(i, p, rp)
                g = Gspan(graph=p)
                dfs = g.minDFS().human_encode()
                pair_idx[p] = dfs
                pair_idx[rp] = dfs

                if dfs not in motif.db:
                    motif.add_with_len(dfs, i)
                    count[i] += 1
                    cumul += 1
                    # parents(pair_idx, motif, p)
                    # print(f'unique {p}  {rp}  {dfs}')
                # else:
                # print(f'\tnot unique {p}  {rp}  {dfs}')

        print(f'level:{i}\tpossible:{possible}\tassymetric:{assymetric}\tunique:{count[i]}\tcumulative'
              f':{cumul}')

    # store as pickled binary file
    ptest = open('pickletest.pkl', 'wb')
    motif.pickle(ptest)
    ptest.close()

    pt = MotifDB.unpickle('pickletest.pkl')

    exit(0)
