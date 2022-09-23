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

    exit(0)
