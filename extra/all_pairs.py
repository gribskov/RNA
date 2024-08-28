"""=================================================================================================
Try to more efficiently generate all structures based on pair representation

Michael Gribskov     18 January 2022
================================================================================================="""
import time


def all_pairs(n):
    """---------------------------------------------------------------------------------------------
    Generate all connected graphs with n stems as a pair list
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
            newpairs = pairs[:] + [left, right]
            # newavail = avail[:].remove(right)
            newavail = avail[:]
            newavail.remove(right)
            if newavail:
                # check if the graph is/will be disconnected, the graph is disconnected all
                # positions less than the current length have been used,
                # i.e. avail[-1] == next available position
                if len(newpairs) == newavail[-1]:
                    # print(f'\tdisconnected {newpairs})')
                    continue

                else:
                    # connected, push on stack
                    stack.append([newpairs, newavail])
            else:
                # reached the desired size, don't push
                yield newpairs


def reverse(pairs):
    """-----------------------------------------------------------------------------------------
    Reverse the pairlist and put in canonical order

    note: faster to pass l
          passing r or creating with list comprehension does not

    :param pairs: list, stem pairs
    :return: list
    -----------------------------------------------------------------------------------------"""
    r = pairs[:]
    l = len(r)
    pos = 0
    for i in sorted(range(1, l, 2), key=lambda x: pairs[x], reverse=True):
        r[pos], r[pos + 1] = l - pairs[i] - 1, l - pairs[i - 1] - 1
        pos += 2

    return r


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    total = 0
    n = 10
    start_time = time.process_time()
    # stage = start_time
    for i in range(1, n):
        count = 0
        stage = time.process_time()
        # r = [0 for _ in range(i*2)]
        # l = i * 2
        for p in all_pairs(i):
            # rp = reverse(p, l)
            rp = reverse(p)
            if rp > p:
                # print(i, p, rp, 'R')
                pass
            else:
                # print(i, p, rp)
                count += 1
            # count += 1

        now = time.process_time()
        total += count
        try:
            persec = count / (now - stage)
        except ZeroDivisionError:
            persec = 0

        print(f'\t{i:2d}\t{count:10d}\t{total:10d}\t{int(persec):10d}\t{now - stage:.2e}\t{now - start_time:.2e}')

    exit(0)
