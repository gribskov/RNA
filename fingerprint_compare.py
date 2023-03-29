"""=================================================================================================
compare old fingerprint to new fingerprint

Michael Gribskov     20 September 2022
================================================================================================="""
import glob
import os
import sys
import argparse
from lxml import etree
from fingerprint import Fingerprint


def process_command_line():
    """---------------------------------------------------------------------------------------------

    :return: namespace      command line arguments
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Compare sets of XIOS fingerprints',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
    )
    cl.add_argument('-n', '--new',
                    help='New format fingerprint (default%(default))',
                    action='append', nargs='*',
                    default='*.fpt')
    cl.add_argument('-h', '--hex',
                    help='Old hexadecimal encoded format fingerprint (default%(default))',
                    action='append', nargs='*',
                    default='*.xpt')
    cl.add_argument('-o', '--outputdir',
                    help='Directory for fingerprint output (default=%(default)s)',
                    default='./')

    args = cl.parse_args()

    return args

def read_new_fpt(target_list):
    """---------------------------------------------------------------------------------------------
    Read a set of fingerprints. Each glob from the filename is put in a separate set of motif names
    in the string format

    in the new format (2022-2023, yaml file0 such as
    - fingerprint:
    - information:
          Date: '2022-05-06 09:27:07'
          File: ./fpt/rnasep_a4.Pseudoanabaena_sp.PCC6903.w4.d5.fpt
          Motif database: ../RNA/data/2to7stem.mdb.pkl
          RNA structure: ./xios/rnasep_a4.Pseudoanabaena_sp.PCC6903.w4.d5.xios
      - total: 50001
      - nmotif: 909
      - motif:
          0i1.0i2.0i3.0o4.: 981
          0i1.1o2.2j0.0i3.0i4.: 877
          ...

    :param target_list: list    string with paths to sets of files to compare
    :return: dict               sets of motifs; set is a list of strings such as 0i1.0i2.0i3.0o4
    ---------------------------------------------------------------------------------------------"""
    motif_set = {}
    for target in target_list:
        fptfilelist = fastafiles = glob.glob(target)
        set = {}

        n = 0
        for fptfile in fptfilelist:
            n += 1
            id = os.path.basename(fptfile)

            # read new fingerprint as YAML
            fpt = Fingerprint()
            fpt.readYAML(fptfile)
            set[id] = fpt.motif

        motif_set[target] = set
        return motif_set

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
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
    opt = process_command_line()

    print(f'new format fingerprints: {}')


    print(f'fingerprint directory: {newdir}')
    print(f'target: {olddir}')

    fptfilelist = fastafiles = glob.glob(target)
    # print(fptfilelist)
    # print()

    n = 0
    comp = {}
    for newfile in fptfilelist:
        oldfile = os.path.basename(newfile)
        id = oldfile.replace('.fpt', '')
        comp[id] = {}
        n += 1
        print(f'{n:4d}\t{id}')
        # print(comp)

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
        new_miss = 0
        new_list = []
        for code in fpt.motif:
            if code not in oldlist:
                new_miss += 1
                new_list.append(code)

        comp[id]['new_miss'] = new_miss
        # comp[id]['new_list'] = new_list

    line = 0
    for fpt in sorted(comp, key=lambda f: ((comp[f]['count_old'] - comp[f]['count_new']), f)):
        line += 1
        # print(fpt, comp[fpt])
        f = comp[fpt]
        print(
            f"{line:4d}\t{f['count_old']}\t{f['count_new']}\t{f['old_miss']}\t{f['new_miss']}\t{fpt}")
        if f['old_list']:
            for motif in f['old_list']:
                print(f'\t\t\t\t\t{motif}')

    exit(0)
