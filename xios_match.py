"""=================================================================================================
compare a set of xios files to a set of reference files and calculate overlap in terms of stems and
bases

Michael Gribskov     23 September 2022
================================================================================================="""
import sys
import os
import glob

from lxml import etree


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

        # use list comprhension to convert fields 3-6 to int and store
        stems.append([int(field[i]) for i in range(3, 7)])

    return stems


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # read reference files
    reffilelist = glob.glob(sys.argv[1])
    refdata = {}
    ref_n = 0
    for ref in reffilelist:
        name = os.path.basename(ref).replace('.xios', '')
        refdata[name] = read_stems_xml(ref)
        ref_n += 1

    print(f'reference topologies: {ref_n}')
    for r in refdata:
        print(f'{r}')
        for m in refdata[r]:
            print(f'\t{m}')

    # for each target file, compare to reference and calculate overlap

    exit(0)
