from topology import Topology, PairRNA
from xios import Gspan


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
        print(f'\tnotunique {s}')
        unique[s] += 1
        unique[r] += 1
        return False

    # if unseen store both forward and reverse, this minimizes the number of reversals
    unique[s] = 1
    unique[r] = 1
    topo.reverse()  # change back to original orientation

    return True


if __name__ == '__main__':

    unique = {}  # index of unique pair graphs
    ucount = [0]
    gparent = {}
    parentdfs = None
    gchild = {}

    level_max = 3
    stack = [[] for _ in range(level_max + 2)]

    # initialize stack with a single stem
    topo = PairRNA(inlist=[0, 0])
    stack[0].append(topo)
    print(topo)

    parent_level = 0
    child_level = parent_level + 1
    while True:

        if stack[parent_level]:
            # if the stack is not empty, pop a parent topology from the stack
            parent = stack[parent_level].pop()
        else:
            # otherwise, increment the level
            parent_level = child_level
            ucount.append(0)
            child_level += 1
            if parent_level > level_max:
                # quit when maximum level is reached
                break
            else:
                # if max is not reached, go head and pop
                parent = stack[parent_level].pop()

        g = Gspan(graph=parent)
        parentdfs = g.minDFS()
        parentstr = str(parentdfs)
        print(f'level: {parent_level}     parent: {parent}     dfs: {parentstr}')

        for left in range(len(parent) * 2 + 1):
            for right in range(left + 1, len(parent) * 2 + 2):

                # add the new stem in all possible locations
                child = parent.duplicate()
                for pair in child.pairs:
                    # add 1 or 2 to the previous stem coordinates to account for
                    # the new stem being added
                    for j in (0, 1):
                        if pair[j] >= left:
                            pair[j] += 1
                        if pair[j] >= right:
                            pair[j] += 1

                child.push_pair([left, right])
                child.reorder()
                print(f'\t{child}', end=' ')

                if not child.connected():
                    # save unconnected graphs to stack (because their children can become connected later)
                    stack[child_level].append(child)
                    continue

                if check_uniqueness(unique, child):
                    stack[child_level].append(child)
                    print(f'\tunique {child}', end=' ')
                    ucount[parent_level] += 1

                    g = Gspan(graph=child)
                    childdfs = g.minDFS()
                    childstr = str(childdfs)
                    print(f'dfs:{childstr}')
                    if parentstr in gparent:
                        # parent is known
                        if childstr not in gparent[parentstr]:
                            gparent[parentstr].append(childstr)
                    else:
                        # unknown parent
                        gparent[parentstr] = [childstr]

                    if childstr in gchild:
                        # known child
                        if parentstr not in gchild[childstr]:
                            gchild[childstr].append(parentstr)
                    else:
                        # unknown child
                        gchild[childstr] = [parentstr]




        # end of loop over adding possible new stems

        print(ucount)

    exit(0)
