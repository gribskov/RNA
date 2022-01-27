from topology import Topology, PairRNA
from xios import Gspan, MotifDB
from githash import get_git_revision_short_hash


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

def parents(pair_idx, motif, pair):
    """---------------------------------------------------------------------------------------------
    Find all parents by deleting one stem from pair, canonicalizing and looking up in pair_idx

    :param pair_idx: dict, key is PairRNA string, value is minDFS string
    :param motif: MotifDB
    :param pair: PairRNA, current structure

    TODO: add minlevel instead of just assuming starting at level 3

    :return:
    ---------------------------------------------------------------------------------------------"""
    parentlist = []
    if len(pair.pairs) < 6:
        # need at least three stems to have an edge in the parent
        return parentlist

    new = pair.duplicate()
    for i in range(0, len(pair.pairs), 2):
        new.pairs = pair.pairs[0:i]+pair.pairs[i+2:]
        new.canonical()
        if new.connected() :
            p = pair_idx[str(new)]
            # print(f'{i}  {pair}: {new}  {new}  {p}')
            if p not in parentlist:
                parentlist.append(p)

    return parentlist

# ===================================================================================================
# main
# ===================================================================================================
if __name__ == '__main__':

    pair_idx = {}  # index of unique pair graphs

    level_max = 8
    level_min = 2
    motif = MotifDB()

    count = [0] * level_min
    cumul = 0

    for i in range(level_min, level_max + 1):
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
                pair_idx[str(p)] = dfs
                pair_idx[str(rp)] = dfs

                if dfs not in motif.db:
                    motif.add_with_len(dfs, i)
                    count[i] += 1
                    cumul += 1
                    motif.parent[dfs] = parents(pair_idx, motif, p)
                    # print(f'{p}     {dfs}     {motif.parent[dfs]}')

        print(f'level:{i}\tpossible:{possible}\tassymetric:{assymetric}\tunique:{count[i]}\tcumulative'
              f':{cumul}')

    # store as pickled binary file
    motif.setdate()
    motif.setname(f'{level_min} to {level_max} stem motifs')
    githash = get_git_revision_short_hash()
    motif.information['source'] = f'make_xios_db:v{githash}'
    filename = f'{level_min}to{level_max}stem.mdb.pkl'
    motif.information['file'] = filename
    pkl = open(filename, 'wb')
    motif.pickle(pkl)
    pkl.close()

    print(motif.toJSON())

    exit(0)
