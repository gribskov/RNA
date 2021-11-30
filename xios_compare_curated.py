"""=================================================================================================
compare xios file to xios of a curated file and evaluate overlap

================================================================================================="""
import sys
from topology import Topology

def summary( xios, title ):
    '''---------------------------------------------------------------------------------------------
    Summarize the stems in the XIOS structure
    :param xios: Topology object
    : param title: string, describes file in summary
    :return: True
    ---------------------------------------------------------------------------------------------'''
    print( f'\n{title} XIOS file of {xios.sequence_id}')
    print( f'\tstems: {len(xios.stem_list)}')
    for s in xios.stem_list:
        print( f'\t{s["name"]}\t{s["left_begin"]}\t{s["left_end"]}\t{s["right_begin"]}\t{s["right_end"]}')

    return True

if __name__ == '__main__':

    refname = sys.argv[1]
    subjectname = sys.argv[2]

    # read reference xios (gold standard)
    ref = Topology()
    ref.XIOSread(refname)
    summary(ref, 'Curated')

    # read subject xios (test structure to compare)
    subject = Topology()
    subject.XIOSread(subjectname)
    summary(subject, 'Test')

    exit(0)


