"""=================================================================================================
some xios files fail during reading with Topology:XIOSread()

Michael Gribskov     11 January 2022
================================================================================================="""
from topology import Topology

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    file = 'data/xiosfiles/g2_b.Escherichia_coli.ATBDi1.w1.d1.xios'
    test = Topology()
    test.XIOSread(file)
    test.XIOSwrite()


    exit(0)
