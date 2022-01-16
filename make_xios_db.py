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
    gparent = {}
    gchild = {}

    level_max = 3
    stack = [[] for _ in range(level_max + 1)]
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
            child_level += 1
            if parent_level > level_max:
                # quit when maximum level is reached
                break
            else:
                # if max is not reached, go head and pop
                parent = stack[parent_level].pop()

        print(f'level: {parent_level}     parent: {parent}')

        for left in range(len(parent) * 2 + 1):
            for right in range(left + 1, len(parent) * 2 + 2):
                child = parent.duplicate()
                # print(f'left:{left}     right:{right}')
                # pair_new = [left, right]

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
                if check_uniqueness(unique, child):
                    # g = Gspan(child.pairs)
                    serial = child.to_SerialRNA()
                    stack[child_level].append(child)

                print(f'\t{child}')

    exit(0)
