"""=================================================================================================
compare xios file to xios of a curated file and evaluate overlap

================================================================================================="""
import sys
from topology import Topology


def summary(xios, title):
    """---------------------------------------------------------------------------------------------
    Summarize the stems in the XIOS structure
    :param xios: Topology object
    : param title: string, describes file in summary
    :return: True
    ---------------------------------------------------------------------------------------------"""
    print(f'\n{title} XIOS file of {xios.sequence_id}')
    print(f'\tstems: {len(xios.stem_list)}')
    for s in xios.stem_list:
        print(
            f'\t{s["name"]}\t{s["left_begin"]}\t{s["left_end"]}\t{s["right_begin"]}\t{s["right_end"]}')

    return True


def check_left(x1, x2, p1, p2):
    """---------------------------------------------------------------------------------------------
    check if the left half stem of structure 1 overlaps a left half stem in structure 2, starting
    at stem p1 and p2 in structure 1 and 2, respectively.
    Stems in structures are sorted by left_begin

    :param x1: Topology object, structure 1
    :param x2: Topology object, structure 2
    :param p1: stem of interest in structure 1
    :param p2: starting stem in structure 2
    :return: list of int, stems in structure 2 whos left half-stems overlap p1's
    ---------------------------------------------------------------------------------------------"""
    overlap_set = []
    stem1 = x1.stem_list[p1]
    pos2 = p2
    while True:
        stem2 = x2.stem_list[pos2]
        if stem2['left_begin'] > stem1['left_end']:
            # end condition for while loop, stem2 is to the right of stem1
            break

        if stem2['left_end'] < stem1['left_begin']:
            # stem2 is to the left of stem1
            pos2 += 1
            continue

        # if stem 2 is not to the left or right it overlaps
        overlap_set.append(pos2)
        # print('\tleft: stem1:{} overlaps stem2:{}'.format(p1, pos2))
        pos2 += 1

    return overlap_set


def check_right(x1, x2, p1, overlap):
    """---------------------------------------------------------------------------------------------
    overlap contains a list of stems in structure 2 that overlap stem p1 in structure 1 on the
    left. Check each one to see if it overlaps on the right

    :param x1: Topology object, structure 1
    :param x2:Topology object, structure 2
    :param p1: int, stem of interest in structure 1
    :param overlap: list of int, overlapping left half stems for stem p1 in structure 1
    :return: list of int, stems in structure 2 that overlap with both left and right halfstems of p1
    ---------------------------------------------------------------------------------------------"""
    both = []
    stem1 = x1.stem_list[p1]
    for p2 in overlap:
        stem2 = x2.stem_list[p2]

        if stem2['right_begin'] > stem1['right_end'] or \
                stem2['right_end'] < stem1['right_begin']:
            # no overlap
            continue

        # print('\tright: stem1:{} overlaps stem2:{}'.format(p1, p2))
        both.append(p2)

    return both


def stat_stem(match1, nstem_2):
    """---------------------------------------------------------------------------------------------
    for the matched stems, calculate the stem level precision and recall.  The index of the matched
    list is the number of the stem in structure 1, the value is a list of overlapping stems in 
    structure 2
    
    :param matched: list of list, structure 2 stems that max the index stem in structure 1
    :param nstem_2: number of stems in structure 2
    :return: 
    ---------------------------------------------------------------------------------------------"""
    # set up a match list for structure 2
    match2 = [[] for i in range(nstem_2)]

    # count true positives and false negatives in structure one while populating the second
    # match list
    tp1 = 0
    fn1 = 0
    for s1 in range(len(match1)):
        if match1[s1]:
            tp1 += 1
            for s2 in match1[s1]:
                match2[s2].append(s1)

        else:
            fn1 += 1

    tp2 = 0
    fn2 = 0
    for s2 in range(len(match2)):
        if match2[s2]:
            tp2 += 1
        else:
            fn2 += 1
                
    recall_1 = tp1 / len(match1)
    precision_1 = tp2 / len(match2)
    recall_2 = tp2 / len(match2)
    precision_2 = tp1 / len(match1)
    print(f'structure 1: recall={recall_1:.3f}\tprecision={precision_1:.3f}')
    print(f'structure 2: recall={recall_2:.3f}\tprecision={precision_2:.3f}')
    print(f'jaccard={(tp1+tp2)/(len(match1)+len(match2)):.3f}')
    return True


def overlap2(x1, x2):
    """---------------------------------------------------------------------------------------------
    compare the stems in structure 1 and structure 2 and identify the overlapping stems
    :param x1: Topology object, structure 1
    :param x2: Topology object, structure 2
    :return:
    ---------------------------------------------------------------------------------------------"""
    p2 = 0
    left = []
    both = []
    for p1 in range(len(x1.stem_list)):
        p2 = 0
        left = check_left(x1, x2, p1, p2)
        if left:
            # one or more overlapping left half-stems
            both.append(check_right(x1, x2, p1, left))
        else:
            continue

        # if both:
        #     for p2 in both:
        #         s1 = x1.stem_list[p1]
        #         s2 = x2.stem_list[p2]
        #         print('overlap x1:{} x2:{}  -  {}\t{}\t{}\t{}  -  {}\t{}\t{}\t{}'.format(
        #             p1, p2,
        #             s1['left_begin'], s1['left_end'], s1['right_begin'], s1['right_end'],
        #             s2['left_begin'], s2['left_end'], s2['right_begin'], s2['right_end']
        #             ))

    # print(both)
    stat_stem(both, len(x2.stem_list))

    return True


if __name__ == '__main__':
    refname = sys.argv[1]
    subjectname = sys.argv[2]

    # read reference xios (gold standard)
    ref = Topology()
    ref.XIOSread(refname)
    summary(ref, 'Curated')

    # read subject xios (test structure to compare)
    subject = Topology()
    subject.XIOSread(subjectname)
    summary(subject, 'Test')

    overlap2(ref, subject)

    exit(0)
