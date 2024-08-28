"""=================================================================================================
ct files stochastically sampled from the partition file seem to have the information needed to find
the pseudoknots (based on a sample of 1 RNaseP).  I think a DeBruijn graph like approach can be used
to trace the alternate structures
================================================================================================="""

def ct_read(ctfile):
    """---------------------------------------------------------------------------------------------
    return the next ct structure in the file, example of format below.

      360  rnasep_a2.Neisseria_meningitidis.fa
    1 C       0    2    0    1
    2 G       1    3    0    2
    3 G       2    4  354    3
    4 G       3    5  353    4
    5 A       4    6  352    5
    6 C       5    7  351    6
    from https://rna.urmc.rochester.edu/Text/File_Formats.html#CT
    line 0:
        Start of first line: number of bases in the sequence
        End of first line: title of the structure
    Each of the following lines, one/base
        Base number: index n
        Base (A, C, G, T, U, X)
        Index n-1
        Index n+1
        Number of the base to which n is paired. No pairing is indicated by 0 (zero).
        Natural numbering. RNAstructure ignores the actual value given in natural numbering, so it is easiest to repeat n here.
    :param ct: file     ct file open for reading
    :return: dict       ct information
    ---------------------------------------------------------------------------------------------"""
    ct = {}
    line_n = 0

    # line 0
    line = ctfile.readline()
    token = line.rstrip().split()
    ct['len'] = token[0]
    ct['seqname'] = token[1]

    # other lines, stop when seqlen is reached
    while True:
        line = ctfile.readline()
        token = line.rstrip().split()


if __name__ == '__main__':
    # read the ct file and create a list for each base with its possible base pairing partners and
    # their counts
    ctfilename = 'data/partition.stochastic.ct'
    ct = open(ctfilename, 'r')


