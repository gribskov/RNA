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
        #the left location in the next pair is always the next lowest location
        left = avail.pop()

        for right in avail:
            # all other locations are available for the right position, push all on the stack
            newpairs = pairs[:] + [left, right]
            # newavail = avail[:].remove(right)
            newavail = avail[:]
            newavail.remove(right)
            if newavail:
                next = 1
                con = True
                for i in range(len(newpairs) - 1):
                    if newpairs[i] == next:
                        print(f'\tdisconnected {newpairs})')
                        con = False
                        break
                    else:
                        next = max(next, newpairs[i] + 1)

                if con:
                    stack.append([newpairs, newavail])
            else:
                print(newpairs)

    exit(0)
