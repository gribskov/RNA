"""=================================================================================================
compar old fingerprint to new fingerprint

Michael Gribskov     20 September 2022
================================================================================================="""
from lxml import etree

def decodedfs( hexstr ):
    """---------------------------------------------------------------------------------------------
    decode the compressed hexadecimal dfs (perl version) to the current python string version

    :param hexstr: string
    :return: string
    ---------------------------------------------------------------------------------------------"""
    xios = ['i','j','o']
    v1mask = 224
    v2mask = 28
    emask = 3
    dfs = ''

    i = 0
    while i < len(hexstr):
        hexword = int(hexstr[i:i+2],16)
        v1 = (hexword & v1mask) >> 5
        v2 = (hexword & v2mask) >> 2
        edge = hexword  & emask
        dfs += f'{v1}{xios[edge]}{v2}.'
        i += 2

    return dfs
# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # get old new fingerprint name amd derive old fingerprint name

    olddir = 'data/fpt/oldfpt/'
    newdir = '/.'

    newfile = '5S_b.Thermus_thermophilus.w4.d5.fpt'
    oldfile = newfile.replace('w4.d5.fpt', 'xios.xpt')

    # read old fingerprint as xml
    xptfile = open(olddir + oldfile, 'r')
    xpt = etree.parse(xptfile)
    print(etree.tostring(xpt))

    codelist = xpt.xpath('//encoded_dfs')
    dfslist = []
    for code in codelist:
        print(f'{code.text}')
        dfs = decodedfs(code.text)
        # dfs = decode('ff')
        print(f'dfs: {dfs}')
        dfslist.append(dfs)

    print(f'{oldfile}\t{len(dfslist)} motifs')

    exit(0)
