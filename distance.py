"""=====================================================================================================================
distance.py
 structural distances

Michael Gribskov 4/16/2026
====================================================================================================================="""
import sys
import pandas as pd


class Distance:
    def __init__(self, filename=' fingerprint.distance'):
        self.filename = filename
        self.df = None

    def read(self):
        """-----------------------------------------------------------------------------------------
        Read distances from self.fh.  Distances are stored in a pandas dataframe in self.df

        :return:
        -----------------------------------------------------------------------------------------"""
        self.df = pd.read_csv(self.filename, delim_whitespace=True)

        return self.df


# ======================================================================================================================
# main
# ======================================================================================================================
if __name__ == '__main__':
    distance_file = 'data/fingerprint.distance'
    dist = Distance(distance_file)
    dist.read()

    print(f'{len(dist.df)} distances read from {distance_file}')

    exit(0)
