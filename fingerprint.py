class Fingerprint(list):
    """=============================================================================================
    A finger print is a spectru of fixed size motifs identified in a structure.  The motifs are
    referenced to a motif dictionary to enable simple identification of parents.

    10 July 2019     Michael Gribskov
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        super().__init__(self)
        self.information = { 'structure_id':'',
                             'motif_db_id':'',
                             'date': ''
                            }
        self.length = 0

# ##################################################################################################
# Testing
# ##################################################################################################
if __name__ == '__main__':

    finger = Fingerprint()
    finger.append('dummy')
    finger.append('dummier')
    finger.append('dummiest')
    print(finger[0])
    print(finger)


    exit(0)
