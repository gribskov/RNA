####################################################################################################
# ct2xios
#
# convert ct file to XIOS formatted XML file
####################################################################################################
import sys
from datetime import datetime
from topology import RNAstructure

# ===================================================================================================
# main
# ===================================================================================================
if __name__ == '__main__':
    ctfile = sys.argv[1]
    sys.stderr.write('CT file: {}'.format(ctfile))

    rna = RNAstructure()
    rna.CTRead(ctfile)
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    rna.comment.append('creation_date {}'.format(now))
    rna.comment.append('input_file {}'.format(ctfile))
    if rna.energy:
        rna.comment.append('input_format unifold')
    else:
        rna.comment.append('input_format probknot')
    rna.stemlist_from_pairs(unpaired=2)
    rna.adjacency_from_stemlist()
    rna.edgelist_from_adjacency(include="ijo", whole=False)

    rna.XIOSwrite(sys.stdout)

    exit(0)
