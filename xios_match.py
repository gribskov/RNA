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
                    help='file glob specifying path to reference files (%(default)s)',
                    default='/depot/mgribsko/rna/curated/curated_Jiajie_Huang_20160220/xios_graph/*.xios')
    cl.add_argument('target_files', nargs='?',
                    help='file glob specifying path to target files (%(default)s)',
                    default='xiosfiles/*.xios')

    return cl.parse_args()


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
    for f in filelist:
        xios = Topology(xios=f)
        stems = []
        for s in xios.stem_list:
            stems.append([s.lbegin, s.lend, s.rbegin, s.rend])
        key = os.path.basename(f)
        key = key.replace('.xios', '')
        data[key] = stems

    return data


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
            # search the target stems for the one that best matches the references stem
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


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
    print(f'match_xios.py {runstart}\n')

    opt = get_options()
    print(f'reference files: {opt.reference_files}', end='')
    refdata = read_xios_stems(opt.reference_files)
    print(f'\t{len(refdata)} files read')

    print(f'target files: {opt.target_files}', end='')
    targetdata = read_xios_stems(opt.target_files)
    print(f'\t{len(targetdata)} files read')
    print()

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

        result = stem_compare(refdata[base], targetdata[target])
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
        print(f'{cond}\t',
              f'\t{result["stem_precision"]:.3f}',
              f'\t{result["stem_recall"]:.3f}',
              f'\t{result["stem_f1"]:.3f}',
              f'\t{result["base_precision"]:.3f}',
              f'\t{result["base_recall"]:.3f}',
              f'\t{result["base_f1"]:.3f}'
              )

    exit(0)
