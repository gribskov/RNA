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
import datetime
import yaml


def process_command_line():
    """---------------------------------------------------------------------------------------------

    :return: namespace      command line arguments
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Compare sets of XIOS fingerprints',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=60)
        )
    cl.add_argument('-n', '--new',
                    help='New format fingerprint (default=%(default)s)',
                    nargs='*', default='*.fpt')
    cl.add_argument('-e', '--encode',
                    help='Old hexadecimal encoded format fingerprint (default=%(default)s)',
                    nargs='*', default='*.xpt')
    cl.add_argument('-o', '--outputdir',
                    help='Directory for fingerprint output (default=%(default)s)',
                    default='./')

    args = cl.parse_args()

    return args


def read_fingerprints(opt):
    """---------------------------------------------------------------------------------------------
    Read old format fingerprint (type='encode) and new fingerprint format (type='new') that match
    the file globs in opt.new and opt.encode. fingerprints are stored in separate sets, one set for
    each file glob on the command line

    :param opt: namespace       command line option namespace from process_command_line()
    :return: list               contents of list are sets of fingerprints, each set is a list of dicts
    ---------------------------------------------------------------------------------------------"""
    fptset = []
    for fpttype in ('encode', 'new'):
    # for fpttype in ('new', 'encode'):
        source = opt.new
        reader = read_new_fpt
        if fpttype == 'encode':
            source = opt.encode
            reader = read_encode_fpt
        for target in source:
            # source is a list of file globs identifying different sets of fingerprints
            this_set = []
            files = glob.glob(target)
            fptset.append(reader(files))

    return fptset

def read_encode_fpt(target_list):
    """---------------------------------------------------------------------------------------------
    read a list of fingerprints in the old (<2022) format, also known as .xpt format, .xpt format
    is XML such as

    <?xml version="1.0"?>
    <XIOS_fingerprint>
        <query>
            <query_id>/home/huang147/reactor/holder20160220/5S_a.Halobacterium_salinarum.xios</query_id>
            <query_vertex>5</query_vertex>
            <query_edge>6</query_edge>
        </query>

        <fingerprint>
            <type>random</type>
            <iteration>500000</iteration>
            <program>fingerprint_random.pl v1.1.4.8</program>
            <time_elapsed>317.450853</time_elapsed>
        </fingerprint>
        <database>
            <database_id>/home/huang147/Motif_fingerprint/2_to_7_stems_topologies.storable</database_id>
        </database>

        <motif_list>
            <motif_n>1</motif_n>
            <motif>
                <id>5_276</id>
                <count>500000</count>
                <first_observed>1</first_observed>
                <encoded_dfs>0428410c7081</encoded_dfs>
                <mapping>1</mapping>
            </motif>
            ...

    :param target_list: list    filenames matching --encode fileglob
    :return: list               contents of list are sets of fingerprints, each set is a list of dicts
    ---------------------------------------------------------------------------------------------"""
    fpt_set = {}
    for target in target_list:
        id = os.path.basename(target)
        prefix = name_prefix(id, 4)
        fpt = Fingerprint()
        xptfile = open(target, 'r')
        xpt = etree.parse(xptfile)

        # read the information fields
        info = etree_to_dict(xpt)
        fpt.information = info

        # read and transform the motifs
        for m in xpt.xpath('//motif_list/motif'):
            code = m.find('encoded_dfs')
            motif = decodedfs(code.text)
            count = m.find('count')
            fpt.motif[motif] = int(count.text)

        fpt_set[prefix] = {'target':target, 'fpt':fpt}

    return fpt_set

def etree_to_dict(xpt):
    """---------------------------------------------------------------------------------------------
    read the <query> <fingerprint> and <database> sections of the fingerprint and return as a dict
    {   'query":{   'query_id':
                    'query_vertex':
                    'query_edge':   }
        'fingerprint':{  'type':
                        'iteration':
                        'program':
                        'time_elapsed:  }
        'database_id':
    }

    :param xpt: etree element       etree parsed xml
    :return: dict                   described above
    ---------------------------------------------------------------------------------------------"""
    d = {}
    # query information
    q = xpt.xpath('//query')
    if q:
        d['query'] = {}
        for tag in ('query_id','query_vertex', 'query_edge'):
            info = xpt.xpath('//query/{}'.format(tag))[0].text
            d['query'][tag] = info

    f = xpt.xpath('//fingerprint')
    if f:
        d['fingerprint'] = {}
        for tag in ('type','iteration', 'program', 'time_elapsed'):
            info = xpt.xpath('//fingerprint/{}'.format(tag))[0].text
            d['fingerprint'][tag] = info

    db = xpt.xpath('//database_id')
    if db:
        d['database_id'] = db[0].text

    return d


def read_new_fpt(target_list):
    """---------------------------------------------------------------------------------------------
    Read a set of fingerprints. fingerprints for all files in target list are returned in a dict
    keyed with their shortened names (see name_prefix())

    new fingerprint format (2022-2023), yaml file such as
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
    :return: dict               all fingerprints in the input target list
    ---------------------------------------------------------------------------------------------"""
    fpt_set = {}
    for target in target_list:
        id = os.path.basename(target)
        prefix = name_prefix(id, 4)

        # read new fingerprint as YAML
        fpt = Fingerprint()
        fpt.readYAML(target)
        fpt_set[prefix] = {'target': target, 'fpt': fpt}

    return fpt_set


def name_prefix(name, n):
    """---------------------------------------------------------------------------------------------
    return a shortened name with the first n tokens, with n=4
    rnasep_a4.Pseudoanabaena_sp.PCC6903.w4.d5.fpt becomes
    rnasep_a4_Pseudoanabaena_sp

    :param name: str    name to shorten
    :param n: int       number of token to include in prefix
    :return: str        shortened name
    ---------------------------------------------------------------------------------------------"""
    begin = 0
    tokens = []
    for i in range(len(name)):
        if name[i] in '._':
            tokens.append(name[begin:i])
            begin = i + 1

    return '_'.join(tokens[:n])


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
    test = name_prefix('d', 4)
    print(f'fingerprint_compare.py fingerprints: {runstart}')
    # xptlist = glob.glob(opt.encode)
    fptlist = read_fingerprints(opt)
    # for target in opt.new:
    #     fptlist.append({'type':'new', 'file':glob.glob(target)})
    # for target in opt.encode:
    #     fptlist.append({'type':'encode', 'file':glob.glob(target)})

    # read

    # fpt = read_new_fpt(fptlist)

    n = 0
    comp = {}
    for newfile in fptlist:
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
