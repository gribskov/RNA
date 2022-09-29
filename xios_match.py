"""=================================================================================================
compare a set of xios files to a set of reference files and calculate overlap in terms of stems and
bases

Michael Gribskov     23 September 2022
================================================================================================="""
import sys
import os
import glob
import argparse
import datetime

from lxml import etree

from topology import Topology


def get_options():
    """---------------------------------------------------------------------------------------------
    Get options from command line using argparse

    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Calculate precision/recall for sets of reference and target xios topologies',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('reference_files', nargs='?',
                    help='path to reference file directory (%(default)s)',
                    default='/depot/mgribsko/rna/curated/curated_Jiajie_Huang_20160220/xios_graph/')
    cl.add_argument('target_files', nargs='?',
                    help='path to target file directory (%(default)s)',
                    default='xiosfiles/')

    opts = cl.parse_args()

    # make sure directories end in '/'
    if not opts.reference_files.endswith('/'):
        opts.reference_files = opts.reference_files + '/'
    if not opts.target_files.endswith('/'):
        opts.target_files = opts.target_files + '/'

    return opts


def read_xios_stems(fileglob):
    """---------------------------------------------------------------------------------------------
    Read stem lists in xios XML format from a fileglob specifying a set of files, i.e., *.xios.
    Returns a dictionary with the base filename as the key and a list of [lbegin, lend, rbegin, rend]
    as values. the suffix '.xios' is removed from the key
    
    :param fileglob:string  specifies a path to a set of xios files
    :return:dictionary  
    ---------------------------------------------------------------------------------------------"""
    filelist = glob.glob(fileglob)
    if not filelist:
        sys.stderr.write(f'No files found in ({fileglob})')
        exit(2)

    data = {}
    seqlen = {}
    for f in filelist:
        xios = Topology(xios=f)
        stems = []
        for s in xios.stem_list:
            stems.append([s.lbegin, s.lend, s.rbegin, s.rend])
        key = os.path.basename(f)
        key = key.replace('.xios', '')
        data[key] = stems
        seqlen[key] = len(xios.sequence)

    return [data, seqlen]


def stats(overlap):
    """---------------------------------------------------------------------------------------------
    return precision, recall and F1
    :param overlap: list [lbegin, lend, rbegin, rend]
    :return: list   precision, recall, F1
    ---------------------------------------------------------------------------------------------"""
    recall = overlap['both'] / overlap['rlen']
    precision = overlap['both'] / overlap['tlen']
    f1 = 0.5 * (precision + recall)

    return precision, recall, f1


def stem_overlap(rstem, tstem):
    """---------------------------------------------------------------------------------------------
    tests to see if both half stems overlap,
    :param rstem: list  [lbegin, lend, rbegin, rend]
    :param tstem: list  [lbegin, lend, rbegin, rend]
    :return:
    ---------------------------------------------------------------------------------------------"""
    result = True
    overlap = {'both': 0, 'rlen': 0, 'tlen': 0, 'maxlen': 0}
    for i in [0, 2]:
        minbegin = min(rstem[i], tstem[i])
        minend = min(rstem[i + 1], tstem[i + 1])
        maxbegin = max(rstem[i], tstem[i])
        maxend = max(rstem[i + 1], tstem[i + 1])
        common = minend - maxbegin + 1

        if common < 0:
            result = False
            break
        else:
            overlap['both'] += minend - maxbegin + 1
            overlap['rlen'] += rstem[i + 1] - rstem[i] + 1
            overlap['tlen'] += tstem[i + 1] - tstem[i] + 1
            overlap['maxlen'] += maxend - minbegin + 1

    if result:
        return overlap
    else:
        return []


def stem_compare(refstemlist, targetstemlist):
    """---------------------------------------------------------------------------------------------
    for each reference stem find the biggest overlap in the targets. If multiple stems overlap, but
    are not the best, they contribute to the denominator of the precision since only the best is
    considered correct.

    :param refstemlist:
    :param targetstemlist:
    :return: list
    ---------------------------------------------------------------------------------------------"""
    total = {'both': 0, 'rlen': 0, 'tlen': 0, 'maxlen': 0}

    stem_n = 0
    match_n = 0
    for rstem in refstemlist:
        # print(rstem)
        stem_n += 1
        tlen_total = 0
        f1_best = 0
        overlap_best = {}
        for tstem in targetstemlist:
            # search the target stems for the one that best matches the reference stem
            overlap = stem_overlap(rstem, tstem)
            if overlap:
                precision, recall, f1 = stats(overlap)
                tlen_total += overlap['tlen']

                if f1 > f1_best:
                    f1_best = f1
                    overlap_best = overlap

                # print(f'\toverlap {tstem}: [{recall:.3f},{precision:.3f}]\t{overlap}')

        if overlap_best:
            # will be unset if there is no overlap of rstem to any tstem
            match_n += 1
            # print(f'\tbest: {overlap_best}')
            total['both'] += overlap_best['both']
            total['rlen'] += overlap_best['rlen']
            # total['tlen'] += overlap_best['tlen']
            total['tlen'] += tlen_total
            total['maxlen'] += overlap_best['maxlen']
        else:
            # if not overlap, add tlen
            total['tlen'] += tstem[1] + tstem[3] - tstem[0] - tstem[2] + 2

    bprecision, brecall, bf1 = stats(total)
    sprecision = match_n / len(targetstemlist)
    srecall = match_n / stem_n
    sf1 = 0.5 * (sprecision + srecall)
    # print(f'\nreference: {stem_n}\ttarget: {len(targetstemlist)}\tmatched:{match_n}')
    # print(f'recall: {recall:.3f}\tprecision: {precision:.3f}\tf1:{f1:.3f}\ttotal: {total}')

    return {'stem_precision': sprecision, 'stem_recall': srecall, 'stem_f1': sf1,
            'base_precision': bprecision, 'base_recall': brecall, 'base_f1': bf1}


def base_compare(refstemlist, targetstemlist, seqlen):
    """---------------------------------------------------------------------------------------------
    calculate base precision and recall by constructing an image of the sequence indicating whether
    each base is found in the reference, predicted in the target, and found in an overlapping
    segment.  Use bits to indicate reference (1), target(2), overlap(4).  Possible values are
    1   only in reference
    2   only in target
    3   reference and target but not in overlap
    7   reference and target and in overlap

    :param refstemlist: list    elements are [lbegin, lend, rbegin, rend]
    :param targetstemlist: list    elements are [lbegin, lend, rbegin, rend]
    :param seqlen: int  length of sequence
    :return: dict   stem and base precision, recall, and f1
    ---------------------------------------------------------------------------------------------"""
    map = [0 for _ in range(seqlen+2)]

    r_overlap = 0
    t_overlap = [0 for _ in range(len(targetstemlist))]
    for rstem in refstemlist:
        # print(rstem)
        map_bases(map, rstem, 1)
        firstpass = True
        firstoverlap = True
        t = 0
        for tstem in targetstemlist:

            # search the target stems for stems that overlap each reference stem
            if firstpass:
                map_bases(map, tstem, 2)

            if tstem[0] > rstem[1]:
                # stems are ordered by lbegin
                break

            overlap = [max(rstem[0], tstem[0]), min(rstem[1], tstem[1]),
                       max(rstem[2], tstem[2]), min(rstem[3], tstem[3])]
            if overlap[0] <= overlap[1] and overlap[2] <= overlap[3]:
                if firstoverlap:
                    r_overlap += 1
                    firstoverlap = False

                t_overlap[t] = 1
                map_bases(map, overlap, 4)

            t += 1
        firstpass = False

    s_rec = r_overlap / len(refstemlist)
    s_pre = sum(t_overlap) / len(targetstemlist)
    s_f1 = 0.5 * (s_rec + s_pre)
    b_pre, b_rec, bf1 = map_stats(map)

    return {'stem_precision': s_pre, 'stem_recall': s_rec, 'stem_f1': s_f1,
            'base_precision': b_pre, 'base_recall': b_rec, 'base_f1': bf1}


def map_bases(map, stem, key):
    """---------------------------------------------------------------------------------------------
    mark the bases in stem by OR with key, expect 1=reference, 2=target, 4=overlap
    :param map:list bit image of sequence with one int per base position
    :param stem:list [lbegin, lend, rbegin, rend]
    :param key:int  a binary key for the source of the stem
    :return: None   map is modified
    ---------------------------------------------------------------------------------------------"""
    for pos in range(stem[0], stem[1] + 1):
        map[pos] = map[pos] | key
    for pos in range(stem[2], stem[3] + 1):
        try:
            map[pos] = map[pos] | key
        except IndexError:
            print(f'index error pos={pos} len={len(map)}')

    return


def map_stats(map):
    """---------------------------------------------------------------------------------------------
    in the final map, compare to the keys and count the number of each type.
    :param map: list bit image of sequence with one int per base position
    :return: float, float, float    precision, recall, f1
    ---------------------------------------------------------------------------------------------"""
    reflen = tarlen = olen = 0
    for pos in map:
        if pos & 1:
            reflen += 1
        if pos & 2:
            tarlen += 1
        if pos & 4:
            olen += 1

    recall = olen / reflen
    precision = olen / tarlen
    f1 = 0.5 * (recall + precision)

    return precision, recall, f1


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
    print(f'match_xios.py {runstart}\n')

    dummy = {}
    opt = get_options()
    print(f'reference files: {opt.reference_files + "*.xios"}', end='')
    refdata, dummy = read_xios_stems(opt.reference_files + '*.xios')
    print(f'\t{len(refdata)} files read')

    print(f'target files: {opt.target_files + "*.xios"}', end='')
    targetdata, seqlen = read_xios_stems(opt.target_files + '*.xios')
    print(f'\t{len(targetdata)} files read')

    # for each target file, compare to reference and calculate overlap
    columns = ['stem_precision', 'stem_recall', 'stem_f1',
               'base_precision', 'base_recall', 'base_f1']
    print(f'# target\tcondition\tsprecision\tsrecall\tsf1\tbprecision\tbrecall\tbf1')
    condition_average = {}
    condition_n = {}

    base_old = ''
    for target in targetdata:
        # get the base name that corresponds to the reference
        suffixpos = target.index('.w')  # make sure no target file names contain .w, or this will fail
        base = target[:suffixpos]
        condition = target[suffixpos + 1:]

        if base_old != base:
            base_old = base
            # print()

        # stemresult = stem_compare(refdata[base], targetdata[target])
        result = base_compare(refdata[base], targetdata[target], seqlen[target])
        print(f'{base}\t{condition}\t',
              f'\t{result["stem_precision"]:.3f}',
              f'\t{result["stem_recall"]:.3f}',
              f'\t{result["stem_f1"]:.3f}',
              f'\t{result["base_precision"]:.3f}',
              f'\t{result["base_recall"]:.3f}',
              f'\t{result["base_f1"]:.3f}'
              )

        if condition in condition_average:
            condition_n[condition] += 1
            for tag in columns:
                condition_average[condition][tag] += result[tag]
        else:
            # new condition
            condition_n[condition] = 1
            condition_average[condition] = {t: result[t] for t in columns}

    # end of loop over targets

    print(f'\n# averages by condition')
    print(f'# condition\tsprecision\tsrecall\tsf1\tbprecision\tbrecall\tbf1')
    for cond in sorted(condition_n.keys(), key=lambda x: condition_average[x]["stem_f1"]):
        for col in columns:
            condition_average[cond][col] /= condition_n[condition]

        result = condition_average[cond]
        print(f'{cond:10s}\t',
              f'\t{result["stem_precision"]:.3f}',
              f'\t{result["stem_recall"]:.3f}',
              f'\t{result["stem_f1"]:.3f}',
              f'\t{result["base_precision"]:.3f}',
              f'\t{result["base_recall"]:.3f}',
              f'\t{result["base_f1"]:.3f}'
              )

    exit(0)
