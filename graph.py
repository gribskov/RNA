class Graph:
    """=============================================================================================
    RNA graph class
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Graph constructor
        -----------------------------------------------------------------------------------------"""
        self.structure = ''


def possible(s):
    """
    return a list of the possible next stems.  Two kinds of stems can be added
    1) stems with use zero can be added in numerical order
    2) stems that have use == 1 can be added, but only if the the count of open stems would be > 0

    :param s: stem usage list
    :return: list of possible elements to add
    """
    possible = []
    candidate = []
    first = True
    open = 0
    all = 0
    for i in range(len(s)):
        all += 2 - s[i]
        if s[i] == 0 and first:
            first = False
            possible.append(i)
        elif s[i] == 1:
            candidate.append(i)
            open += 1

    if open > 1 or all == 1:
        possible += candidate

    return possible


def enumerate(n):
    """---------------------------------------------------------------------------------------------
    enumerate all graphs with n stems
    :param n: number of stems
    :return: array of graph lists
    ---------------------------------------------------------------------------------------------"""
    graphs = []
    s = [0 for i in range(n)]
    stack = []
    struc = []
    stack.append((0, s[:], struc))

    while stack:
        n, s, struc = stack.pop()
        struc.append(n)
        s[n] += 1

        p = possible(s)
        if p:
            for n in p:
                stack.append((n, s[:], struc[:]))
        else:
            graphs.append(struc)

    return graphs

# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    total = 0
    for size in range(1, 8):
        g = enumerate(size)
        total += len(g)
        print(size, len(g), total)

graphs = enumerate(2)
for g in graphs:
    print(g)

exit(0)
