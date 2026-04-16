"""=====================================================================================================================
distance_stats.py

calculates statistics on structure distances

Michael Gribskov 4/16/2026
====================================================================================================================="""
import glob
import sys
import os
from collections import defaultdict
from fingerprint import Fingerprint, FingerprintSet
import numpy as np
import pandas as pd


def get_name_group(namestr):
    """-----------------------------------------------------------------------------------------------------------------
    get the structure name and group from the xios file name

    :param namestr: string     xios structure file, from fingerprint yaml
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    name = os.path.basename(namestr)
    name.replace('.xios', '')
    group = name
    if '_' in group:
        group = group.split('_')[0]
    if '.' in group:
        group = group.split('.')[0]
    print(f'{nfpt}\t{name}\t class: {group}')

    return name, group


def entropy(data):
    """-----------------------------------------------------------------------------------------------------------------

    :param data:
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    column_totals = df.sum()


    return


def plottotal(data, size):
    # print(sorted(d, reverse=True))
    import matplotlib.pyplot as plt
    # sort the data
    sorted_data = np.sort(data)[::-1]

    # x positions
    x = np.arange(len(sorted_data))

    # plot
    plt.figure()
    plt.bar(x, sorted_data)
    plt.gca().xaxis.set_visible(False)

    plt.xlabel("Sorted index")
    plt.ylabel("Value")
    plt.title(f"Sorted data bar chart n={size}")
    plt.show()

    return


# ======================================================================================================================
# main
# ======================================================================================================================
if __name__ == '__main__':
    fptglob = '../data/fpt/*.xios.out'

    fpt_list = glob.glob(fptglob)
    print(f'fingerprints: {fptglob}')

    motifs = defaultdict(list)
    fptidx = []
    groupidx = []
    nfpt = 0
    records = []
    for this_fpt in fpt_list:
        fpt = Fingerprint()
        nfpt += 1
        fpt.readYAML(this_fpt)
        name, group = get_name_group(fpt.information['RNA structure'])

        for f in fpt.motif:  # set() avoids double-counting within sample
            records.append({
                "sample" : name,
                "group"  : group,
                "feature": f
            })

    df = pd.DataFrame(records)

    total_counts = df.groupby("feature")["sample"].nunique()
    group_counts = (
        df.groupby(["feature", "group"])["sample"]
        .nunique()
        .unstack(fill_value=0)
    )

    result = group_counts.copy()
    result["total"] = total_counts

    # plottotal(result['total'], nfpt)
    entropy(result)

    exit(0)
