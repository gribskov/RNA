def find_tranche(data):
    """---------------------------------------------------------------------------------------------
    Find data groups that have the same value, referred to here as a tranche

    :param data: list       int or float proxy statistic values
    :return: list           positions of beginning of each tranche
    ---------------------------------------------------------------------------------------------"""
    # tranche = [0]
    # value = data[0]
    # for i in range(1, len(data)):
    #     if data[i] == value:
    #         continue
    #
    #     else:
    #         tranche.append(i)
    #         value = data[i]
    tranche = []
    value = data[0]
    inside = False
    for i in range(1, len(data)):
        if data[i] == value:
            # next value is the same (a tranche)
            if inside:
                continue
            else:
                tranche.append(i - 1)
                inside = True
        else:
            # next value is different (end open tranches)
            if inside:
                inside = False
                tranche.append(i)

        value = data[i]
    if inside:
        tranche.append(len(data))

    return tranche

def tranche2(data):
    tranche = []
    value = data[0]
    dir = label[0]
    p0 = n0 = p = n = 0
    i = 1
    while i < len(data):
        # i is the next value (incremented at end of loop)

        # add the next range of values
        while value == label[i]:
            # when i exit, i will be the first point with a new value
            dir = None
            p += label[i]
            n += not label[i]
            i += 1

        if dir == label[i]:
            # next value is the same continue in current direction
            p1 += label[i]
            n1 += not label[i]
        else:
            # direction changed - calculate area
            a = a2(n, p, p0)
            p0 += p
            n = 0

        # update direction, value, and i

        dir = label[i]
        value = value[i]
        i += 1


def a2(n, p, p0):
    pass

    return

def area(p0, p1, pn, n0, n1, nn):
    """--------------------------------------------------------------------------------------------
    calculate trapezoidal area
    :param p0: int      positive begin
    :param p1: int      positive end
    :param pn: int      n positive
    :param n0: int      negative begin
    :param n1: int      negative begin
    :param nn: int      n negative
    :return: float      area
    --------------------------------------------------------------------------------------------"""
    if n1 = n0:
        return 0.0

    a = 0.5 * (p0 + p1) / pn
    a *= (n1 - n0) / nn

    return a




def roc(data, label, tranche):
    """---------------------------------------------------------------------------------------------

    :param data:
    :param label:
    :return:
    ---------------------------------------------------------------------------------------------"""
    pos_n = label.count(True)
    neg_n = label.count(False)

    current = data[0]
    dir = label[0]
    i = 0
    t = tranche.pop(0)
    p0 = n0 = p1 = n1 = 0
    calc = False
    while i < len(data):
        if i == t:
            # area of current chunk
            a = area(p0, p1, pos_n, n0, n1, neg_n)
            if label[i]:
                p1 = p1 +1
            else:
                n1 = n1 + 1
            p0 = p1
            n0 = n1
            # starting to process tranche
            end = tranche.pop(0)
            while i < end:
                p1 += label[i]
                n1 += not label[i]
                i += 1
            calc = True
            t = tranche.pop(0)
        else:
            if label[i] == dir:
                # continuing in current direction, don't do anything
                p1 += label[i]
                n1 += not label[i]
            else:
                # need to calculate
                calc = True

        if calc:
            # calculate an area
            a = area(p0, p1, pos_n, n0, n1, neg_n)


        i += 1

    return


data = [9, 8, 7, 7, 6, 6, 5, 4, 4, 2, 2]
label = [True, True, True, False, True, True, False, False, True, False, False]

print(data)
tranche = find_tranche(data)
print(tranche)
roc(data, label, tranche)
