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
        # print(f'\tnotunique {s}')
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
    ucount = [0,1]
    gparent = {}
    parentdfs = None
    gchild = {}

    level_max = 4
    # stack = [[] for _ in range(level_max + 2)]

    # initialize stack with a single stem
    topo = PairRNA(inlist=[0, 0])
    stack = [topo]
    ucount = [0]
    # print(topo)
    # parent_level = 0
    # child_level = parent_level + 1
    parent_old = 0
    while True:

        # if parent_level >= level_max:
        #     # if the stack is not empty, pop a parent topology from the stack
        #     break

        try:
            parent = stack.pop()
            parent_level = len(parent)
            child_level = parent_level + 1

        except IndexError:
            print('got here')
            break

        g = Gspan(graph=parent)
        parentdfs = g.minDFS()
        parentstr = str(parentdfs)
        if parent_level > parent_old:
            print(f'level: {parent_level}     parent: {parent}     dfs: {parentstr}')
            parent_old = parent_level
            ucount.append(0)

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
                print(f'\tchecking {child}')

                # if not child.connected():
                #     # save unconnected graphs to stack (because their children can become connected later)
                #     if len(child) < level_max:
                #         stack.append(child)
                #     continue

                if check_uniqueness(unique, child):
                    if len(child) < level_max:
                        stack.append(child)
                    print(f'\t\tunique')
                    ucount[parent_level] += 1

                    if not child.connected():
                        # do not run gspan on disconnected graphs
                        continue

                    print('\t\tconnected')
                    g = Gspan(graph=child)
                    childdfs = g.minDFS()
                    childstr = str(childdfs)
                    print(f'\t\tdfs:{childstr}')
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
    cc = [0]*100
    for cs in gchild:
        for letter in ('[',']',',',"'",' '):
            cs.replace(letter,'')
        l = len(cs)
        print(l, l/3)
        cc[len(cs)] += 1

    print(cc)

    exit(0)
