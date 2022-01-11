"""=================================================================================================
compare xios file to xios of a curated file and evaluate overlap

================================================================================================="""
import sys
import glob
import os.path
from topology import Topology


def summary(xios, title):
    """---------------------------------------------------------------------------------------------
    Summarize the stems in the XIOS structure
    :param xios: Topology object
    : param title: string, describes file in summary
    :return: True
    ---------------------------------------------------------------------------------------------"""
    print(f'\n{title} XIOS file {xios.sequence_id}')
    print(f'\tstems: {len(xios.stem_list)}')
    for s in xios.stem_list:
        print(f'\t{s.name}\t{s.lbegin}\t{s.lend}\t{s.rbegin}\t{s.rend}')

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
        if pos2 >= len(x2.stem_list):
            break

        stem2 = x2.stem_list[pos2]
        if stem2.lbegin > stem1.lend:
            # end condition for while loop, stem2 is to the right of stem1
            break

        if stem2.lend < stem1.lbegin:
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

        if stem2.rbegin > stem1.rend or \
                stem2.rend < stem1.rbegin:
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
    jaccard = (tp1 + tp2) / (len(match1) + len(match2))

    return {'precision': precision_1, 'recall': recall_1, 'jaccard': jaccard}


def overlap(x1, x2):
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

    return both


def parse(filename, refdir):
    """---------------------------------------------------------------------------------------------
    parse the reference name, p, d, and w parameters from the filename
    expects a name like rnasep_m.M_maripaludus.w1.d3.xios

    :param filename: string
    :param refdir: string, path to directory with reference files
    :return: dict, keys = ('name', 'd', 'w' )
    ---------------------------------------------------------------------------------------------"""
    base = os.path.basename(filename)
    field = base.split('.')
    suffix = field.pop()
    d = field.pop().lstrip('d')
    w = field.pop().lstrip('w')
    field += [suffix]

    reference = refdir + '.'.join(field)

    return {'name': reference, 'd': d, 'w': w}


# ===================================================================================================
# main program
# ===================================================================================================

if __name__ == '__main__':
    # refname = sys.argv[1]
    # subjectglob = sys.argv[2]

    xiosdir = sys.argv[1]
    refdir = sys.argv[2]
    if not refdir.endswith('/'):
        refdir += '/'

    # for each xios file, determine the reference name from the xios filename and read
    # for a name like rnasep_m.M_jannaschii.w1.d0.xios
    # the reference would be rnasep_m.M_jannaschii.xios

    current_ref = None
    overall = {}

    for testfile in glob.glob(xiosdir + '/*.xios'):
        # print(testfile)
        parsed = parse(testfile, refdir)

        if os.path.exists(parsed['name']):
            # read reference xios (gold standard)
            if parsed['name'] != current_ref:
                ref = Topology()
                ref.XIOSread(parsed['name'])
        else:
            # skip if there is no reference file
            continue

        ref.sequence_id = parsed['name']
        current_ref = os.path.basename(parsed['name'])
        # split on whatever comes first, . or _
        family = current_ref.split('.')[0]
        family = f'{family.split("_")[0]}.w{parsed["w"]}.d{parsed["d"]}'
        all = f'all.w{parsed["w"]}.d{parsed["d"]}'
        if family not in overall:
            overall[family] = {'precision': 0, 'recall': 0, 'jaccard': 0, 'n': 0}
        if all not in overall:
            overall[all] = {'precision': 0, 'recall': 0, 'jaccard': 0, 'n': 0}
        # summary(ref, 'Curated')

        # read subject xios (test structures to compare)
        subject = Topology()
        subject.XIOSread(testfile)
        # summary(subject, 'Test')

        both = overlap(ref, subject)
        stat = stat_stem(both, len(subject.stem_list))
        overall[family]['precision'] += stat['precision']
        overall[family]['recall'] += stat['recall']
        overall[family]['jaccard'] += stat['jaccard']
        overall[family]['n'] += 1
        overall[all]['precision'] += stat['precision']
        overall[all]['recall'] += stat['recall']
        overall[all]['jaccard'] += stat['jaccard']
        overall[all]['n'] += 1
        # print('{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{}\t{}'.format(
        #     family, parsed['d'], parsed['w'], stat['precision'],
        #     stat['recall'], stat['jaccard'],
        #     len(ref.stem_list), len(subject.stem_list), current_ref))

    for k in sorted(overall):
        (id, w, d) = k.split('.')
        w = int(w[1])
        d = int(d[1])
        n = overall[k]['n']
        print('{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}'.format(
            id, w, d, n,
            overall[k]['precision'] / n, overall[k]['recall'] / n, overall[k]['jaccard'] / n))

    exit(0)
