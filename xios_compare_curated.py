"""=================================================================================================
compare xios file to xios of a curated file and evaluate overlap
================================================================================================="""
import sys
from topology import Topology

if __name__ == '__main__':

    refname = sys.argv[1]
    subjectname = sys.argv[2]

    # read reference xios (gold standard)
    ref = Topology()
    ref.XIOSread(refname)

    # read subject xios
    subject = Topology()
    subject.XIOSread(subjectname)

    exit(0)


