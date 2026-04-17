"""=====================================================================================================================
distance_stats.py

calculates statistics on structure distances

Michael Gribskov 4/16/2026
====================================================================================================================="""
import glob
import math
import sys
import os
#from cgi import print_form
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


def entropy(data, prior):
    """-----------------------------------------------------------------------------------------------------------------

    :param prior: series      fraction of motifs in each group: expected distribution of a random motif
    :param data:dataframe     motif counts, columns are sample groups, rows are motifs
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    # print(data.info())
    # print(data.head())
    # print(total)
    groups = data.columns.tolist()
    groups.remove('all')

    h0 = 0.0
    for group, value in prior.items():
        if group != 'all':
            h0 += -prior[group] * math.log2(prior[group])
    # print(f'background entropy: {h0: .3f}')

    f = data
    f = data[groups].div(data['all'], axis=0)
    h = -(f * np.log2(f)).sum(axis=1) - h0

    print(h.head())

    return h.sort_values(ascending=False)

    def plottotal(data, size):
        """-----------------------------------------------------------------------------------------------------------------

        :param data:
        :param size:
        :return:
        -----------------------------------------------------------------------------------------------------------------"""
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
    #fptglob = '../data/fpt/*.xios.out'
    fptglob = '*.fpt'

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

    group_counts = (
        df.groupby(["feature", "group"])["sample"]
        .nunique()
        .unstack(fill_value=0)
    )

    result = group_counts.copy()
    result['all'] = df.groupby("feature")["sample"].nunique()
    # print(result.head())

    total = result.sum(numeric_only=True)
    prior = total / total['all']
    plus = result +  prior
    print(result.head())
    # plottotal(result['total'], nfpt)
    h = entropy(plus, prior)
    result['entropy'] = h
    print('\ncounts with entropy')
    pd.options.display.max_rows = 2000
    pd.options.display.max_columns = 20
    pd.options.display.max_colwidth = 100
    pd.options.display.width = 200
    print(result.sort_values(by='entropy'))

    exit(0)
