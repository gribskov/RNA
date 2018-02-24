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


def toPairs(g):
    """"--------------------------------------------------------------------------------------------
    convert a graph in array format to pair format.  returns a list of tuples with the begin/end
    position of each stem
    :param g: graph in array format
    :retuen: array of tuples
    ---------------------------------------------------------------------------------------------"""
    pairs = [[] for i in range(int(len(g) / 2))]

    for i in range(len(g)):
        pairs[g[i]].append(i)

    return pairs


def fromPairs(p):
    """---------------------------------------------------------------------------------------------
    convert a graph in pair format to array format
    :param p: graph in pair format
    :return g: graph in list format
    ---------------------------------------------------------------------------------------------"""
    g = [0 for i in range(len(p) * 2)]

    stem = 0
    for pair in p:
        g[pair[0]] = stem
        g[pair[1]] = stem
        stem += 1

    return g


def reversePairs(p):
    """---------------------------------------------------------------------------------------------
    reverse the positions of the stems by converting to maxpos-pos and resorting to be in order of
    first coordinate
    :param p: graph in pair format
    :return: graph in pair format (reversed)
    ---------------------------------------------------------------------------------------------"""
    m = len(p) * 2 - 1
    for i in range(len(p)):
        p[i][0], p[i][1] = m - p[i][1], m - p[i][0]


    p.sort(key=lambda k: k[0])

    return p


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    total = 0
    for size in range(1, 8):
        g = enumerate(size)
        total += len(g)
        print(size, len(g), total)

graphs = enumerate(3)
for g in graphs:
    print('\ngraph:', g)
    pairs = toPairs(g)
    print('pairs:', pairs, end='\t=>\t')
    print(fromPairs(pairs))
    rev = reversePairs(pairs)
    print('reversed:', rev)

exit(0)
