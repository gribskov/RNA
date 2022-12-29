def ROC(data, label):
    """---------------------------------------------------------------------------------------------
    Calculate the receiver operating characteristic (ROC) curve and the area under the curve (AUC).
    This is an exact calculation taking into account identical values.

    Data values and labels (True, False) are provided as two lists. Use the sort option if the data
    is not sorted by value.

    :param data: list of int/float  data values
    :param label: list of bool      data labels
    :param sort: bool               sort the data and labels
    :param sortdir: string          sort direction high/low only the first char is checked, case insensitive
    :return: list of float/ float   points for drawing curve, AUC
    ---------------------------------------------------------------------------------------------"""
    # get numbers of positives and negatives
    pn = nn = 0
    for tf in label:
        pn += tf
        nn += not tf

    dir = label[0]
    value = data[0]
    area = 0
    p0 = n0 = p = n = 0

    p += dir
    n += not dir

    i = 1
    points = [[0, 0]]
    while i < len(data) - 1:

        # look forward, is value the same
        firstpt = True
        while data[i] == value:
            # read in tranche, a tranche is a set of points with the same data value
            # stopping i = first different value
            p += label[i]
            n += not label[i]
            if firstpt:
                # draw to the starting point of the tranche
                points.append([n0 + n, p0 + p])
                firstpt = False

            dir = None  # set direction == None
            if i + 1 < len(data):
                i += 1
            else:
                break

        while label[i] == dir and data[i] != value:
            # will not run after reading a tranche of same value points
            p += label[i]
            n += not label[i]
            value = data[i]
            if i + 1 < len(data):
                i += 1
            else:
                break

            if data[i] == value:
                # back out the first point in a subsequent tranche
                i -= 1
                p -= label[i]
                n -= not label[i]
                value = data[i]

        # at this point, i points to the value that begins the next block and the count buffers
        # are filled with the previous block.

        if points[-1] != [n0 + n, p0 + p]:
            points.append([n0 + n, p0 + p])

        # calculate area

        a = (p0 + 0.5 * p) * n
        area += a
        # print(f'i:{i}\tp0:{p0}\tp:{p}\tn:{n}\ta:{a}\t{area}')

        # update count buffers
        p0 += p
        n0 += n
        p = n = 0
        dir = label[i]
        value = data[i]
        # print(f'\tp0:{p0}\tp:{p}\tn:{n}\tdir:{dir}\tvalue:{value}')

    # end of loop over data points

    if points[-1] != [n0 + n, p0 + p]:
        points.append([n0 + n, p0 + p])

    if nn:
        # if number of negatives > 0
        area /= pn * nn

    for pt in points:
        pt[0] /= nn
        pt[1] /= pn

    return points, area


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    data = [9, 8, 6, 6, 6, 6, 5, 4, 4, 2, 2]
    label = [True, True, True, False, True, True, False, False, True, False, False]
    # data = [1, 1, 2, 2, 3, 3.5, 4, 4, 5, 5]
    # label = [False, False, False, False, False, True, True, True, True, True, ]

    for i in range(len(data)):
        print(f'{i}  {data[i]}  {label[i]}')

    # tranche = find_tranche(data)
    points, auc = ROC(data, label)
    print(f'\nfinal area:{auc:.4f}')
    print(points)

    exit(0)
