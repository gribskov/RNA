"""=================================================================================================


Michael Gribskov     15 December 2022
================================================================================================="""
import sys


def ROC2(roc_file, distance, score, npos, nneg):
    """---------------------------------------------------------------------------------------------
    calculate Receiver Operating Characteristic and area under curve (AUC)

    :param roc_file:    string, filename for ROC output
    :param distance:    list of dict, distance information
    :param npos:        int, number of positive comparisons
    :param nneg:        int, number of negative comparisons
    :return:
    ---------------------------------------------------------------------------------------------"""

    nstep = 1.0 / nneg
    pstep = 1.0 / npos
    pbeg = pend = 0
    nbeg = nend = 0
    dirup = True
    value = 0
    auc = area = 0.0
    for point in sorted(distance, key=lambda x: x[score], reverse=True):
        print(f'pend={pend}/{npos}\tnend={nend}/{nneg}\tvalue={value}\t{dirup}')
        if point[score] == value:
            # items with same value cannot be sorted so they are a single step
            if point['ispos']:
                pend += 1
            else:
                nend += 1
            continue

        else:
            # value has changed but you don't need to calculate until both pos
            # and negative values change
            if not dirup and not point['ispos']:
                # this step is negative, following previous negative, just increment neg
                nend += 1
            if dirup and point['ispos']:
                # positive step following previous positives, i.e., straight up
                pend += 1
                pbeg += 1
            else:
                # area of trapezoid
                area = (pbeg + pend) * pstep * (nend - nbeg) * nstep / 2.0
                auc += area
                roc_file.write(
                    f'{pend * pstep:.6f}\t{nend * nstep:.6f}\t{point[score]}\t{area:.4g}\t{auc:.6f}{newline}')
                # print(f'score:{point[score]}\t area:auc:')

                if point['ispos']:
                    pend += 1
                    dirup = True
                else:
                    dirup = False
                    # nend += 1

                pbeg = pend
                nbeg = nend

            value = point[score]
            print(f'    pend={pend}/{npos}\tnend={nend}/{nneg}\tvalue={value}\t{dirup}')

    roc_file.write(f'{pend * pstep:.4f}\t{(nend) * nstep:.4f}\t{point[score]}\t{area:.4f}\t{auc:.4f}{newline}')
    print(f'pend={pend}/{npos}\tnend={nend}/{nneg}')
    area = (pbeg + pend) * pstep * (nend - nbeg) * nstep / 2.0
    auc += area

    return auc


def roc3(out, distance, score, npos, nneg):
    d = sorted(distance, key=lambda x: x[score], reverse=True)

    dir = d[0]['ispos']
    prev = None
    p = 0
    n = 0
    i = 0
    pbeg = 0
    pend = 0
    nbeg = 0
    nend = 0
    area = 0.0
    trapezoid = False

    for i in range(len(d)):
        print(f'start: {pbeg} {pend}\t{nbeg} {nend}\t{d[i][score]}\t{area}')
        pos = d[i]['ispos']
        pend += pos
        nend += not pos
        if d[i][score] == prev:
            trapezoid = True
            continue

        if nend - nbeg > 0:
            if trapezoid:
                area += 0.5 * (pend + pbeg) / npos * (nend - nbeg) / nneg
                print(f'  trapezoid {pbeg} {pend}\t{nbeg} {nend}\t {area}')
                trapezoid = False
            else:
                area += pend / npos * (nend - nbeg ) / nneg
                print(f'  rectangle {pbeg} {pend}\t{nbeg} {nend}\t {area}')
            pbeg = pend
            nbeg = nend

        prev = d[i][score]
        print(f'  end: {pbeg} {pend}\t{nbeg} {nend}\t {area}\n')

    return

def ranges(distance):
    # ranges end at the end of a negative run
    range = []
    beg = 0
    end = 0
    negative = False
    for i in range(len(d)):
        if d[i]['ispos']:
            if negative:
                range.append((beg,i))
                beg = i + 1
            else:
                negative = True

    return

def clump(distance):
    score = 'score'
    d = sorted(distance, key=lambda x: x[score], reverse=True)
    clump = []

    beg = 0
    # end = 0

    intranche = False
    if d[0][score] == d[1][score]:
        intranche = True

    value = d[0][score]
    sign = d[0]['ispos']
    i = 0
    while i < len(d)-1:
        nextvalue = d[i+1][score]
        nextsign = d[i+1]['ispos']

        if d[i][score] == nextvalue:
            # current score is same as next one
            # if intranche:
            #     end = i
            # else:
            if not intranche:
                # not in a tranche, end region
                clump.append([beg, i-1])
                beg = i
                # end = i + 1
                intranche = True
        else:
            # current score is different from next one
            if intranche:
                clump.append([beg,i])
                beg = i + 1
                intranche = False
            # else:
            #     end += 1

        i += 1

    if clump[-1][1] != len(d)-1:
        clump.append([beg, i])

    print(clump)
    return

def clump2(distance):
    score = 'score'
    d = sorted(distance, key=lambda x: x[score], reverse=True)
    npos = nneg = 0
    different = [0]
    previous = d[0][score]
    for i in range(1,len(d)):
        current = d[i][score]
        if d[i]['ispos']:
            npos += 1
        else:
            nneg += 1
        if  current != previous:
            different.append(i)
        previous = current

    pbeg = nbeg = 0
    np = nn = 0
    i = 0
    while i < len(d):
        stop = different.pop(0)
        while i < stop:
            np += d[i]['ispos']
            nn += not d[i]['ispos']
            i += 1

        if nn == 0:
            pass
        else:
            print(f'calculate {i}')
            if np > 0 and nn > 0:
                a = (2 * pbeg + np)/npos * (nn)/nneg
                np = nn = 0
            elif np > 0:
                pass
            elif nn > 0:
                a = pbeg/npos * nn/nneg
                np = nn = 0
            pbeg = i

    print(different)
    return






# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
newline = '\n'
distance = []
for i in range(0, 5):
    distance.append({'score': 10 - i, 'ispos': True})
for i in range(5, 10):
    distance.append({'score': 10 - i, 'ispos': False})

scores = [10, 10, 8, 8, 8, 5, 4, 3, 2, 1]
pos    = [1,  1, 1, 0, 1, 1, 0, 0, 0, 0]
npos = 0
nneg = 0
for i in range(len(scores)):
    distance[i] = {'score':scores[i], 'ispos': bool(pos[i])}
    if pos[i]:
        npos += 1
    else:
        nneg += 1

# roc3(sys.stdout, distance, 'score', npos, nneg)
clump2(distance)