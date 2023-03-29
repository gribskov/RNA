"""=================================================================================================
compare old fingerprint to new fingerprint

Michael Gribskov     20 September 2022
================================================================================================="""
import glob
import os
import sys

from lxml import etree

from fingerprint import Fingerprint


def decodedfs(hexstr):
    """---------------------------------------------------------------------------------------------
    decode the compressed hexadecimal dfs (perl version) to the current python string version. The
    hexadecimal dfs represents each row of the dfs as

    000 000 000
    v1  v2  edge where the edge values are 00=i 01=j 10=o

    :param hexstr: string
    :return: string
    ---------------------------------------------------------------------------------------------"""
    xios = ['i', 'j', 'o']
    v1mask = 224
    v2mask = 28
    emask = 3
    dfs = ''

    i = 0
    while i < len(hexstr):
        hexword = int(hexstr[i:i + 2], 16)
        v1 = (hexword & v1mask) >> 5
        v2 = (hexword & v2mask) >> 2
        edge = hexword & emask
        dfs += f'{v1}{xios[edge]}{v2}.'
        i += 2

    return dfs


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # get old new fingerprint name amd derive old fingerprint name

    # olddir = 'data/fpt/oldfpt/'
    # newdir = 'data/fpt/'
    olddir = sys.argv[1]
    newdir = sys.argv[2]
    target = newdir + '*.fpt'

    print(f'fingerprint directory: {newdir}')
    print(f'target: {olddir}')

    fptfilelist = fastafiles = glob.glob(target)
    #print(fptfilelist)
    #print()

    n = 0
    comp = {}
    for newfile in fptfilelist:
        oldfile = os.path.basename(newfile)
        id = oldfile.replace('.fpt','')
        comp[id] = {}
        n += 1
        print(f'{n:4d}\t{id}')
        #print(comp)

        # read new fingerprint as YAML
        fpt = Fingerprint()
        fpt.readYAML(newfile)
        comp[id]['count_new'] = len(fpt.motif)

        # read old fingerprint as xml
        oldfile = oldfile.replace('.fpt', '.xios.xpt')
        xptfile = open(olddir + oldfile, 'r')
        xpt = etree.parse(xptfile)
        oldxpt_to_decode = xpt.xpath('//encoded_dfs')

        oldlist = []
        for code in oldxpt_to_decode:
            motif = decodedfs(code.text)
            oldlist.append(motif)

        comp[id]['count_old'] = len(oldlist)

        # see if old motifs exist in new
        old_miss = 0
        old_list = []
        for code in oldlist:
            if code not in fpt.motif:
               old_miss += 1
               old_list.append(code)

        comp[id]['old_miss'] = old_miss
        comp[id]['old_list'] = old_list

        # see if new motifs exist in old
        new_miss =  0
        new_list = []
        for code in fpt.motif:
            if code not in oldlist:
                new_miss += 1
                new_list.append(code)

        comp[id]['new_miss'] = new_miss
        #comp[id]['new_list'] = new_list


    line = 0
    for fpt in sorted(comp, key=lambda f: ((comp[f]['count_old']-comp[f]['count_new']),f)):
        line += 1
        #print(fpt, comp[fpt])
        f = comp[fpt]
        print(f"{line:4d}\t{f['count_old']}\t{f['count_new']}\t{f['old_miss']}\t{f['new_miss']}\t{fpt}")
        if f['old_list']:
            for motif in f['old_list']:
                print(f'\t\t\t\t\t{motif}')

    exit(0)
