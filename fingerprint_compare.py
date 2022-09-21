"""=================================================================================================
compar old fingerprint to new fingerprint

Michael Gribskov     20 September 2022
================================================================================================="""
import glob
import os
from lxml import etree
from fingerprint import Fingerprint

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
    newdir = 'data/fpt/'
    target = newdir+'*.fpt'

    fptfilelist = fastafiles = glob.glob(target)

    for newfile in fptfilelist:
        # newfile = '5S_b.Thermus_thermophilus.w4.d5.fpt'
        oldfile = os.path.basename(newfile)
        oldfile = oldfile.replace('w4.d5.fpt', 'xios.xpt')

        # read new fingerprint as YAML
        fpt = Fingerprint()
        fpt.readYAML(newfile)
        newfile = os.path.basename(newfile)

        # read old fingerprint as xml
        xptfile = open(olddir + oldfile, 'r')
        xpt = etree.parse(xptfile)
        # print(etree.tostring(xpt))

        codelist = xpt.xpath('//encoded_dfs')
        dfslist = []
        for code in codelist:
            # print(f'{code.text}')
            dfs = decodedfs(code.text)
            # dfs = decode('ff')
            # print(f'dfs: {dfs}')
            dfslist.append(dfs)

        print(f'{oldfile}\t{len(dfslist)} motifs')
        print(f'{newfile}\t{len(fpt.motif)} motifs\n')

    exit(0)
