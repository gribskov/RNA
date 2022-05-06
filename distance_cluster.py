"""=================================================================================================
Read the distances calculated by fingerprint_distance. Expected format
fpt1    fpt2    Jaccard Bray-Curtis

First cluster in to connected components, then run UPGMA on each cluster

Michael Gribskov     06 May 2022
================================================================================================="""


def read_distance(filename, maximum, minimum):
    """---------------------------------------------------------------------------------------------
    Read distances from fingerprint_distance.py.  Distances are store in a list of dict
    distance = [ {fpt1, fpt2, jaccard, bray:curtis}, ...]
    
    :param filename: string, name of file with distances
    :return: distance (list of dict), maximum (float), minimum (float)
    ---------------------------------------------------------------------------------------------"""
    pass

    return distance, maximum, minimum


def connected(distance, threshold):
    """---------------------------------------------------------------------------------------------
    Read distances from fingerprint_distance.py

    :param distance: list of dict, keys: fpt1, fpt2, jaccard, bray-curtis
    :return: list of string, IDs of fingerprints in connected components
    ---------------------------------------------------------------------------------------------"""
    pass

    return


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------

distance, max, min = read_distance(distance_file)

components = connected(distance, threshold)

if __name__ == '__main__':
    exit(0)
