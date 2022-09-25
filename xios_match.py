"""=================================================================================================
compare a set of xios files to a set of reference files and calculate overlap in terms of stems and
bases

Michael Gribskov     23 September 2022
================================================================================================="""
import sys
import os
import glob

from lxml import etree

from topology import Topology


def read_stems_xml(refname):
    """---------------------------------------------------------------------------------------------

    :param refname:
    :return:
    ---------------------------------------------------------------------------------------------"""
    reffile = open(refname, 'r')
    xpt = etree.parse(reffile)

    stemlist = xpt.xpath('//stem_list')
    text = stemlist[0].text
    stems = []
    for line in text.split('\n'):
        field = line.split()
        if len(field) < 1:
            # blank line
            continue

        # use list comprehension to convert fields 3-6 to int and store
        stems.append([int(field[i]) for i in range(3, 7)])

    return stems


def read_reference(fileglob):
    """---------------------------------------------------------------------------------------------
    Read the reference topologies in sios XML format and store as a dictionary of
    lists. Each dictionary value is a list of the stems given as a list of
    [left_begin, left_end, right_begin, right_end]

    :param fileglob: string     reference topology files (XML format)
    :return: dictionary
    ---------------------------------------------------------------------------------------------"""
    # read reference files
    reffilelist = glob.glob(fileglob)
    refdata = {}
    for ref in reffilelist:
        name = os.path.basename(ref).replace('.xios', '')
        refdata[name] = read_stems_xml(ref)

    return refdata


def read_target(fileglob):
    """---------------------------------------------------------------------------------------------
    Read target files in xios XML format
    # TODO since reference and target are both in XML, they can be read as xios objects

    :param fileglob:
    :return:
    ---------------------------------------------------------------------------------------------"""
    targetlist = glob.glob(fileglob)
    targetdata = {}
    for target in targetlist:
        xios = Topology(xios=target)
        stems = []
        for s in xios.stem_list:
            stems.append([s.lbegin, s.lend, s.rbegin, s.rend])
        targetdata[target] = stems

    return targetdata


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
        f1_best = 0
        overlap_best = {}
        for tstem in targetstemlist:
            overlap = stem_overlap(rstem, tstem)
            if overlap:
                precision, recall, f1 = stats(overlap)

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
            total['tlen'] += overlap_best['tlen']
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
    refdata = read_reference(sys.argv[1])
    print(f'\treference topologies: {len(refdata)}')

    #     for m in refdata[r]:
    #         print(f'\t{m}')

    # for each target file, compare to reference and calculate overlap
    targetdata = read_target(sys.argv[2])

    columns = ['stem_precision', 'stem_recall', 'stem_f1',
               'base_precision', 'base_recall', 'base_f1']
    condition_average = {}
    condition_n = {}

    base_old = ''
    for target in targetdata:
        # get the base name that corresponds to the reference
        name = os.path.basename(target)
        suffixpos = name.index('.w')  # make sure no names have a .w, or this will fail
        base = name[:suffixpos]
        condition = name[suffixpos+1:len(name)-5]

        if base_old != base:
            base_old = base
            print()
            # print(f'\n{base}')
        # print(f'target:{target}\treference:{base}\texists:{base in refdata}')

        result = stem_compare(refdata[base], targetdata[target])
        print(f'{base}\t{condition}\t',
              f'\tstem precision:{result["stem_precision"]:.3f}',
              f'\tstem recall:{result["stem_recall"]:.3f}',
              f'\tstem f1:{result["stem_f1"]:.3f}',
              f'\tbase precision:{result["base_precision"]:.3f}',
              f'\tbase recall:{result["base_recall"]:.3f}',
              f'\tbase f1:{result["base_f1"]:.3f}'
              )

        if condition in condition_average:
            condition_n[condition] += 1
            for tag in columns:
                condition_average[condition][tag] += result[tag]
        else:
            # new condition
            condition_n[condition] = 1
            condition_average[condition] = {t:result[t] for t in columns}

    # end of loop over targets

    print(f'averages by condition\n')
    for cond in condition_n.keys():
        for col in columns:
            condition_average[cond][col] /= condition_n[condition]

        result = condition_average[cond]
        print(f'{cond}\t{condition}\t',
              f'\tstem precision:{result["stem_precision"]:.3f}',
              f'\tstem recall:{result["stem_recall"]:.3f}',
              f'\tstem f1:{result["stem_f1"]:.3f}',
              f'\tbase precision:{result["base_precision"]:.3f}',
              f'\tbase recall:{result["base_recall"]:.3f}',
              f'\tbase f1:{result["base_f1"]:.3f}'
              )


    exit(0)
