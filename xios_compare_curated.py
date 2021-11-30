"""=================================================================================================
compare xios file to xios of a curated file and evaluate overlap

================================================================================================="""
import sys
from topology import Topology


def summary(xios, title):
    '''---------------------------------------------------------------------------------------------
    Summarize the stems in the XIOS structure
    :param xios: Topology object
    : param title: string, describes file in summary
    :return: True
    ---------------------------------------------------------------------------------------------'''
    print(f'\n{title} XIOS file of {xios.sequence_id}')
    print(f'\tstems: {len(xios.stem_list)}')
    for s in xios.stem_list:
        print(
            f'\t{s["name"]}\t{s["left_begin"]}\t{s["left_end"]}\t{s["right_begin"]}\t{s["right_end"]}')

    return True


def overlap(x1, x2, p1, p2):
    '''---------------------------------------------------------------------------------------------
    return overlap on per stem (True/False) and per base (int) basis
    advance stem position p1 or p2 as needed

    each stem in x1.stemlist and x2.stemlist is a dict of
        {'name', 'center', 'left_begin', 'left_end', 'right_begin', 'right_end', 'left_vienna', 
        'right_vienna' }

    :param x1: Topology object, XIOS structure 1
    :param x2: Topology object, XIOS structure 2
    :param p1: int, current stem position in structure 1
    :param p2: int, current stem position in structure 2
    :return: boolean, int - stem overlaps, number of overlapping bases
    ---------------------------------------------------------------------------------------------'''
    stem = False
    base = 0

    sl_1 = x1.stem_list
    sl_2 = x2.stem_list
    if sl_1[p1]['left_begin'] > sl_2[p2]['left_end']:
        print(f'{p1} > {p2}')
        p2 += 1
        return (p1, p2, stem, base)

    if sl_1[p1]['left_end'] < sl_2[p2]['left_begin']:
        print(f'{p1} < {p2}')
        p1 += 1
        return (p1, p2, stem, base)

    # only reach her if left half stem overlaps, now check the right side

    if sl_1[p1]['right_begin'] > sl_2[p2]['right_end']:
        print(f'{p1} > {p2}')
        p2 += 1
        return (p1, p2, stem, base)

    if sl_1[p1]['right_end'] < sl_2[p2]['right_begin']:
        print(f'{p1} < {p2}')
        p2 += 1
        return (p1, p2, stem, base)

    # both sides overlap
    print('{}|{}\t{}\t{}\t{}\t\t{}|{}\t{}\t{}\t{}'.format(
        p1, sl_1[p1]['left_begin'], sl_1[p1]['left_end'],
        sl_1[p1]['right_begin'], sl_1[p1]['right_end'],
        p2, sl_2[p2]['left_begin'], sl_2[p2]['left_end'],
        sl_2[p2]['right_begin'], sl_2[p2]['right_end']
        ))
    stem = True
    p2 += 1

    return p1, p2, stem, base


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

    p1 = 0
    p2 = 0
    while p1 < len(ref.stem_list):
        p1, p2, stem, base = overlap(ref, subject, p1, p2)

    exit(0)
