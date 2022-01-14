####################################################################################################
# ct2xios
#
# convert ct file to XIOS formatted XML file. Output is STDOUT so redirect to a filename
#
# Usage:
#   ct2xios <ct_file> > <XIOS file>
####################################################################################################
import sys
import argparse
from datetime import datetime
from topology import RNAstructure

def formatter(prog):
    """---------------------------------------------------------------------------------------------
    Set up formatting for help
    :param prog:
    :return: argparse formatter class
    ---------------------------------------------------------------------------------------------"""
    return argparse.HelpFormatter(prog, max_help_position=50, width=120)


def options():
    """---------------------------------------------------------------------------------------------
    Command line options parsed with argparse

    :return: argparse Namespace (similar to a dictionary)
    ---------------------------------------------------------------------------------------------"""
    commandline = argparse.ArgumentParser(
        description='Convert CT file to XIOS file',
        formatter_class=formatter)

    commandline.add_argument('ctfile')
    commandline.add_argument('-d', '--ddG',
                             help='delta deltaG limit for supoptimal structures (%(''default)s) ',
                             type=int,
                             default='5')

    args = commandline.parse_args()
    return args

# ===================================================================================================
# main
# ===================================================================================================
if __name__ == '__main__':

    args = options()
    # sys.stderr.write('CT file: {}\n'.format(args.ctfile))

    rna = RNAstructure()
    rna.CTRead(args.ctfile, ddG=args.ddG)
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    rna.comment.append('creation_date {}'.format(now))
    rna.comment.append('input_file {}'.format(args.ctfile))
    if rna.energy:
        rna.comment.append('input_format unifold')
    else:
        rna.comment.append('input_format probknot')
    rna.adjacency_from_stemlist()
    rna.edgelist_from_adjacency(include="ijo", whole=False)

    rna.XIOSwrite(sys.stdout)

    exit(0)
