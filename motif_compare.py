"""#################################################################################################
Count the number of fingerprints in which each motif is found so that rare or non-discriminatory
motifs can be excluded from the distance calculation

all fingerprints matching fptsuffix are used
#################################################################################################"""
import glob
import os
import sys
from datetime import datetime


class Motif:
    """=============================================================================================
    List of motifs per finger print file
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.motif = {'all': {}}
        self.list = []
        self.n = 0

        return

    def motifread(self, target):
        """-----------------------------------------------------------------------------------------
        Read all the fingerprint files and return a dictionary indexed my fpt name and motif with
        the counts. In addition, an extra fpt, 'all', is returned with the total count of each
        motif.

        motif[fpt name][motifname] = count of motif in fpt

        the values returned in motif are the sampled counts of each motif

        :param target: string       glob to use as target file specification
        :return: dict               motif
        -----------------------------------------------------------------------------------------"""
        motif = self.motif
        for fpt in glob.glob(target):
            # outerloop over fingerprint files
            base = os.path.basename(fpt)
            motif[base] = {}

            fptin = open(fpt, 'r')
            # loop over motifs in fingerprint
            for line in fptin:
                # skip past the header and find the motif section
                if line.find('- motif') > -1:
                    break

            for line in fptin:
                field = line.split(':', maxsplit=1)
                motifname = field[0]
                count = int(field[1].rstrip())

                # motif name and count
                motif[base][motifname] = count
                if motifname in motif['all']:
                    motif['all'][motifname] += count
                else:
                    motif['all'][motifname] = count
                    self.list.append(motifname)
                    self.n += 1

            fptin.close()

        return motif

    def icountgroup(self, groupfxn=None):
        """-----------------------------------------------------------------------------------------
        returns a dictionary with the number of fingerprints in which each motif occurs.
        fxn(fpt_name) transforms the fingerprint name to a group name. an additional group 'all'
        counts the number of occurrences across all groups

        :param groupfxn: function   returns a group name
        :return: dict, dict         key=group, motif value= #of containing fingerprints
                                    key=group, number of fingerprints in the group
        -----------------------------------------------------------------------------------------"""
        if not groupfxn:
            # default group function returns the fingerprint name
            groupfxn = self.nop

        icount = {}  # number of fingerprints that have a motif, by group
        ncount = {'all': 0}  # number of fingerprints in group
        fptlist = self.motif
        for fpt in fptlist:
            if fpt == 'all':
                # skip the all entry
                continue

            group = groupfxn(fpt)
            if group in icount:
                # group already exists
                ncount[group] += 1
            else:
                # first time group has been seen
                icount[group] = {}
                ncount[group] = 1
            ncount['all'] += 1

            for motif in fptlist[fpt]:
                if motif in icount[group]:
                    icount[group][motif] += 1
                else:
                    icount[group][motif] = 1

        # go back and get the total counts
        grouplist = list(icount.keys())
        total = icount['all'] = {}
        for group in grouplist:
            for motif in icount[group]:
                if motif in total:
                    total[motif] += icount[group][motif]
                else:
                    total[motif] = icount[group][motif]

        return icount, ncount

    @staticmethod
    def nop(fptname):
        """-----------------------------------------------------------------------------------------
        icountgroup callback function that does nothing (for use as default)

        :return: string fptname
        -----------------------------------------------------------------------------------------"""
        return fptname


# ===================================================================================================
# end of class Motif
# ===================================================================================================

def process_command_line():
    """---------------------------------------------------------------------------------------------
    read command line options and return as a Namespace object. The Namespace object behaves much as
    a dictionary, and can be converted to a dictionary using the vars() method

    :return: Namespace object
    ---------------------------------------------------------------------------------------------"""
    import argparse

    cl = argparse.ArgumentParser(
        description='Compare and select motifs',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
    )
    cl.add_argument('-p', '--prefix',
                    help='Directory with fingerprint files, can be wildcard (default=%(default)s)',
                    default='../data/fpt/')
    cl.add_argument('-s', '--suffix',
                    help='filename suffix for fingerprint files, can be wildcard  (default=%(default)s',
                    default='*.fpt')
    cl.add_argument('-m', '--motif',
                    help='Distance output file name (default=%(default)s, auto=calculate from fpt)',
                    default='selected.motiflist')
    cl.add_argument('-c', '--cutoff', nargs='*',
                    help='Min Max number of motif occurrences (default=%(default)s)', type=int,
                    default=[])

    args = cl.parse_args()

    if type(args.cutoff) == 'list':
        pass

    for path in ('prefix',):
        thispath = getattr(args, path)
        if not thispath.endswith('/'):
            thispath += '/'
            setattr(args, path, thispath)

    return args


def curated_group(name):
    """---------------------------------------------------------------------------------------------
    icountgroup callback function
    return the name of a group based on the fingerprint name. For the curated files the group is the
    first part of the filename ending at the firs '.' or '_'

    :param name:        fingerprint filename
    :return: string     group name
    ---------------------------------------------------------------------------------------------"""
    field = name.split('.', maxsplit=1)
    field = field[0].split('_')

    return field[0]


def select_by_minmax(icount, cutoff):
    """---------------------------------------------------------------------------------------------
    select a set of motifs that meet a minimum threshold of occurrences (cutoff) across all motifs

    :param icount: dict             count matrix from Motif.icountgroup
    :param cutoff: list of int      counts must >= cutoff[0] and <= cutoff[1]
    :return: dict                   motif names, counts
    ---------------------------------------------------------------------------------------------"""
    selected = {}
    for motif in icount['all']:
        if icount['all'][motif] >= cutoff[0] and icount['all'][motif] <= cutoff[1]:
            selected[motif] = icount['all'][motif]

    return selected


# ===================================================================================================
# main
###=================================================================================================
if __name__ == '__main__':
    newline = '\n'
    opt = process_command_line()
    sys.stderr.write(f'motif_compare: compare and select motifs for use{newline}')
    sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}{newline}{newline}")
    sys.stderr.write(f"target fingerprints: {opt.prefix}{opt.suffix}{newline}")
    sys.stderr.write(f"selected motifs: {opt.motif}{newline}")
    sys.stderr.write(f"cutoff: {opt.cutoff}{newline}")

    try:
        out = open(opt.motif, 'w')
    except OSError:
        sys.stderr.write(f"{newline}Error opening selected motif file ({opt.motif}{newline}")
        exit(1)

    motif = Motif()
    motif.motifread(f"{opt.prefix}{opt.suffix}")

    icount, ncount = motif.icountgroup(curated_group)

    # min and max cutoff
    minmax = []
    if len(opt.cutoff) == 0:
        minmax = [0, ncount['all']]
    elif len(opt.cutoff) == 1:
        minmax = [opt.cutoff[0], ncount['all']]
    else:
        minmax = [opt.cutoff[0], opt.cutoff[1]]

    selected = select_by_minmax(icount, minmax)
    sys.stderr.write(f"{newline}{len(selected)} motifs with counts >= {minmax} written to {opt.motif}{newline}")
    for motif in selected:
        out.write(f'{motif}\t{selected[motif]}{newline}')

    out.close()

    exit(0)
