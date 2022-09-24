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


def stem_overlap(rstem, tstem):
    """---------------------------------------------------------------------------------------------
    tests to see if both half stems overlap,
    :param rstem: list  [lbegin, lend, rbegin, rend]
    :param tstem: list  [lbegin, lend, rbegin, rend]
    :return:
    ---------------------------------------------------------------------------------------------"""
    result = True
    overlap = []
    for i in [0, 2]:
        minbegin = min(rstem[i], tstem[i])
        minend = min(rstem[i + 1], tstem[i + 1])
        maxbegin = max(rstem[i], tstem[i])
        maxend = max(rstem[i + 1], tstem[i + 1])
        overlap.append({'both':   minend - maxbegin + 1,
                        'rlen':   rstem[i + 1] - rstem[i] + 1,
                        'tlen':   tstem[i + 1] - tstem[i] + 1,
                        'maxlen': maxend - minbegin + 1})
        if overlap[-1]['both'] < 0:
            result = False
            break

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
    total = { 'both':0, 'rlen':0, 'tlen':0, 'maxlen':0 }

    for rstem in refstemlist:
        print(rstem)
        f1_best = 0
        overlap_best = []
        for tstem in targetstemlist:
            overlap = stem_overlap(rstem, tstem)
            if overlap:
                recall = (overlap[0]['both'] + overlap[1]['both'])
                recall /= (overlap[0]['rlen'] + overlap[1]['rlen'])
                precision = (overlap[0]['both'] + overlap[1]['both'])
                precision /= (overlap[0]['tlen'] + overlap[1]['tlen'])
                f1 = 0.5 * (precision + recall)

                if f1 > f1_best:
                    f1_best = f1
                    overlap_best = overlap

                print(f'\toverlap {tstem}: [{recall:.3f},{precision:.3f}]\t{overlap}')

        if overlap_best:
            # will be unset if there is no overlap of rstem to any tstem
            print(f'best: {overlap_best}')
            total['both'] += overlap_best['both']
            total['rlen'] += overlap_best['rlen']
            total['tlen'] += overlap_best['tlen']
            total['maxlen'] += overlap_best['maxlen']

    print(f'\ttotal: {total}')

    return


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

    for target in targetdata:
        # get the base name that corresponds to the reference
        name = os.path.basename(target)
        suffixpos = name.index('.w')  # make sure no names have a .w, or this will fail
        base = name[:suffixpos]
        print(f'target:{target}\treference:{base}\texists:{base in refdata}')
        result = stem_compare(refdata[base], targetdata[target])

    exit(0)
