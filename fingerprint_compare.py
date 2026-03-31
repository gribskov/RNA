"""=====================================================================================================================
Compare the motifs in fingerprints one directory to ones with the same name in another directory. This is useful for
comparing fingerprints after changing the calculation

Repurposes old version that compared old fingerprint (XML) to new fingerprint (YAML)

Michael Gribskov     20 September 2022
====================================================================================================================="""
import glob
import os
import sys
import argparse
# from lxml import etree
from fingerprint import Fingerprint
import datetime
# import yaml


def process_command_line():
    """---------------------------------------------------------------------------------------------
    Get command line options
    :return: namespace      command line arguments
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Compare XIOS fingerprints in source to target',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=60)
    )
    cl.add_argument('-s', '--source',
                    help='Source fingerprint glob (default=%(default)s)',
                    nargs='*', default='*.fpt')
    cl.add_argument('-t', '--target',
                    help='Target fingerprint directory (default=%(default)s)',
                    nargs='*', default=[])

    args = cl.parse_args()

    return args


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    daytime = datetime.datetime.now()
    runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')

    opt = process_command_line()
    sys.stderr.write(f'fingerprint_compare.py fingerprints: {runstart}\n')
    sys.stderr.write(f'\tsource fingerprints: {opt.source}\n')
    sys.stderr.write(f'\ttarget fingerprint directory: {opt.target}\n')

    for source_fpt in glob.glob(opt.source[0]):
        sfpt = Fingerprint()
        sfpt.readYAML(source_fpt)

        base = os.path.basename(source_fpt)
        target_fpt = opt.target[0] + base
        tfpt = Fingerprint()
        tfpt.readYAML(target_fpt)

        print(f'comparing {base}')
        present = []
        missing = []
        tmissing = []
        for motif in sfpt.motif:
            if motif in tfpt.motif:
                present.append(motif)
            else:
                missing.append(motif)
        npresent = len(present)
        nmissing = len(missing)
        print(f'\t{len(sfpt.motif)} motifs')
        print(f'\t{npresent} present in {target_fpt}')
        print(f'\t{nmissing} missing in {target_fpt}')

        for motif in tfpt.motif:
            if motif not in sfpt.motif:
                tmissing.append(motif)

        ntmissing = len(tmissing)
        recall = npresent / (npresent + ntmissing)
        precision = npresent / (npresent + ntmissing)
        jaccard = npresent / (npresent + nmissing + ntmissing)
        print(f'\t{ntmissing} missing in {source_fpt}')
        print(f'\trecall: {recall}\tprecision: {precision}\tjaccard: {jaccard}\n')

    exit(0)
