"""=================================================================================================
Read the distances calculated by fingerprint_distance. Expected format
fpt1    fpt2    Jaccard Bray-Curtis

First cluster in to connected components, then run UPGMA on each cluster

Michael Gribskov     06 May 2022
================================================================================================="""
import sys


def read_distance(filename):
    """---------------------------------------------------------------------------------------------
    Read distances from fingerprint_distance.py.  Distances are store in a list of dict
    distance = [ {fpt1, fpt2, jaccard, bray_curtis}, ...]
    
    :param filename: string, name of file with distances
    :return: distance (list of dict), maximum (float), minimum (float)
    ---------------------------------------------------------------------------------------------"""
    try:
        distance = open(filename, 'r')
    except OSError:
        sys.stderr.write(f'distance_cluster - read_distance: cannot open distance file{filename}')

    maximum = {'jaccard': 0, 'bray-curtis': 0}
    minimum = {'jaccard': 1000000, 'bray-curtis': 1000000}

    distance_n = 0
    distance_list = []
    for line in distance:
        if line.startswith('#'):
            continue

        distance_n += 1
        fpt1, fpt2, jaccard, bray_curtis = line.rstrip().split('\t')
        maximum['jaccard'] = max(maximum['jaccard'], float(jaccard))
        minimum['jaccard'] = min(minimum['jaccard'], float(jaccard))
        maximum['bray-curtis'] = max(maximum['bray-curtis'], float(bray_curtis))
        minimum['bray-curtis'] = min(minimum['bray-curtis'], float(bray_curtis))
        distance_list.append({'fpt1':        fpt1,
                              'fpt2':        fpt2,
                              'jaccard':     float(jaccard),
                              'bray-curtis': float(bray_curtis)})

    return distance_list, maximum, minimum


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
distance_file = 'distance.out'
distance, maximum, minimum = read_distance(distance_file)

components = connected(distance, threshold)

if __name__ == '__main__':
    exit(0)
