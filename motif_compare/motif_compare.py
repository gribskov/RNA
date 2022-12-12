"""#################################################################################################
Count the number of fingerprints in which each motif is found so that rare or non discriminatory
motifs can be excluded from the distance calculation

all fingerprints matching fptsuffix are used
#################################################################################################"""
import glob
import os


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


def select_by_min(icount, cutoff):
    """---------------------------------------------------------------------------------------------
    select a set of motifs that meet a minimum threshold of occurences (cutoff) across all motifs

    :param icount: dict     count matrix from Motif.icountgroup
    :param cutoff: int      counts must >= cutoff
    :return: dict           motif names, counts
    ---------------------------------------------------------------------------------------------"""
    selected = {}
    for motif in icount['all']:
        if icount['all'][motif] >= cutoff:
            selected[motif] = icount['all'][motif]

    return selected


"""=================================================================================================
main
================================================================================================="""
if __name__ == '__main__':
    fptglob = 'data/fpt/*'
    fptsuffix = '.fpt'
    cutoff = 8

    motif = Motif()
    motif.motifread(f'{fptglob}{fptsuffix}')

    icount, ncount = motif.icountgroup(curated_group)

    selected = select_by_min(icount, cutoff)
    print(f'{len(selected)} genes with counts >= {cutoff}')
    for motif in selected:
        print(f'{motif}\t\t{selected[motif]}')

    exit(1)
