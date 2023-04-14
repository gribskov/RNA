"""=================================================================================================
ct files stochastically sampled from the partition file seem to have the information needed to find
the pseudoknots (based on a sample of 1 RNaseP).  I think a DeBruijn graph like approach can be used
to trace the alternate structures
================================================================================================="""

if __name__ == '__main__':
    # read the ct file and create a list for each base with its possible base pairing partners and
    # their counts
    ctfilename = 'data/partition.stochastic.ct'
    ct = open(ctfilename, 'r')

