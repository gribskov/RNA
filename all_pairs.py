"""=================================================================================================
Try to more efficiently generate all structures based on pair representation

Michael Gribskov     18 January 2022
================================================================================================="""

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    n = 4
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
                    print(f'\tdisconnected {newpairs})')
                    continue

                else:
                    # connected, push on stack
                    stack.append([newpairs, newavail])
            else:
                # reached the desired size, don't push
                print(newpairs)

    exit(0)
