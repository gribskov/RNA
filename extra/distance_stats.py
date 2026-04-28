"""=====================================================================================================================
distance_stats.py

calculates statistics on structure distances

Michael Gribskov 4/16/2026
====================================================================================================================="""
import glob
import math
import os
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

from fingerprint import Fingerprint


def get_name_group(namestr):
    """-----------------------------------------------------------------------------------------------------------------
    get the structure name and group from the xios file name

    :param namestr: string     xios structure file, from fingerprint yaml
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    name = os.path.basename(namestr)
    name = name.replace('.xios', '')
    group = name
    if '_' in group:
        group = group.split('_')[0]
    if '.' in group:
        group = group.split('.')[0]

    return name, group


def entropy(data, prior):
    """-----------------------------------------------------------------------------------------------------------------
    calculate the entropy of each motif from the frequencies of the motif in each group. Add a pseudocount which is the
    average number of motifs per group summed across all motifs. This corrects for both group size and sequence length.
    the prior distribution is also used to calculate a background entropy which is subtracted

    :param prior: series      fraction of motifs in each group: expected distribution of a random motif
    :param data:dataframe     motif counts, columns are sample groups, rows are motifs
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    # print(data.info())
    # print(data.head())
    # print(total)
    groups = data.columns.tolist()
    groups.remove('all')

    # background entropy
    h0 = 0.0
    for group, value in prior.items():
        if group != 'all':
            h0 += -prior[group] * math.log2(prior[group])

    # entropy for each motif
    f = data
    f = data[groups].div(data['all'], axis=0)
    h = -(f * np.log2(f)).sum(axis=1) - h0
    # print(h.head())

    return h.sort_values(ascending=False)


def binary_entropy(data, prior):
    """-------------------------------------------------------------------------------------------------------------
    two class entropy for distinguishing each group from all others

    :param prior: series      fraction of motifs in each group: expected distribution of a random motif
    :param data:dataframe     motif counts, columns are sample groups, rows are motifs
    :return:
    -------------------------------------------------------------------------------------------------------------"""
    groups = data.columns.tolist()
    groups.remove('all')

    # background entropy
    h0 = pd.Series(np.zeros(len(groups)), index=groups)

    for group, value in prior.items():
        if group != 'all':
            # group vs not-group
            h0[group] += -prior[group] * math.log2(prior[group])
            h0[group] += -(1 - prior[group]) * math.log2(1 - prior[group])

    # entropy for each motif
    f = data
    f = f[groups].div(f['all'], axis=0)

    for group, value in prior.items():
        if group != 'all':
            f[group] = -f[group] * np.log2(f[group]) + (1.0 - f[group]) * np.log2(1.0 - f[group]) - h0[group]
            # h = -(f * np.log2(f)).sum(axis=1) - h0
            # print(h.head())

    col = f.columns
    labels = []
    for c in col:
        labels.append(c + '_be')
    f.columns = labels

    return f


def plottotal(data, size):
    """-------------------------------------------------------------------------------------------------------------
    plot the distribution of the numbr of occurrences of each motif

    :param data:
    :param size:
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
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

def readfpt(fpt_list):
    """-----------------------------------------------------------------------------------------------------------------
    read the kist of fingerprints into a dataframe the raw dataframe is three columns: structure_name, group, and feature
    example:
    16S_b.Clostridium_innocuum     16S     0i1.1i2.2j0.0i3.0i4.0i5.

    :param fpt_list:list    paths to fingerprint files
    :return: dataframe      structure_name, group, feature
    :return: dataframe      count of number of structures in each group
    -----------------------------------------------------------------------------------------------------------------"""
    group = defaultdict(int)
    nfpt = 0
    records = []
    for this_fpt in fpt_list:
        nfpt += 1
        fpt = Fingerprint()
        fpt.readYAML(this_fpt)
        n, g = get_name_group(fpt.information['RNA structure'])
        print(f'\t{nfpt}\t{g}\t{n}')
        group[g] += 1

        for f in fpt.motif:  # set() avoids double-counting within sample
            records.append({
                'sample': n,
                'group': g,
                'feature': f
            })

    raw = pd.DataFrame(records)
    motif = (
        raw.groupby(['feature', 'group'])['sample']
        .nunique()
        .unstack(fill_value=0)
    )

    return raw, group


# ======================================================================================================================
# main
# ======================================================================================================================
if __name__ == '__main__':
    count_min = 5
    count_max = 160
    target_ratio = 0.9
    fptglob = sys.argv[1]

    print(f'fingerprints: {fptglob}')
    print(f'minimum motif count: {count_min}')
    print(f'maximum motif count: {count_max}')
    print(f'target feature fraction: {target_ratio}')

    fpt_list = glob.glob(fptglob)
    print(f'fingerprints: {fptglob}')
    raw, group = readfpt(fpt_list)
    print(f'fingerprints read: {len(raw)}')
    print(f'data groups read: {len(group)}')

    # convert to dataframe with counts of motfis for each group, and filter my maximum and minimum count
    motif = (
        raw.groupby(['feature', 'group'])['sample']
        .nunique()
        .unstack(fill_value=0)
    )
    # group summary: number of samples and motifs per group (before filtering)
    ginfo = pd.DataFrame(data=group, index=['sample_n'])
    ginfo.loc['motif_n'] = motif.sum(axis=0)
    motif_all = ginfo.loc['motif_n'].sum()

    # filter by minimum and maximum count; prior is the probability that a random motif belongs to a group
    motif = motif[(motif.sum(axis=1) > count_min) & (motif.sum(axis=1) < count_max)]
    ginfo.loc['filtered_n'] = motif.sum(axis=0)
    motif_all = ginfo.loc['filtered_n'].sum()
    ginfo.loc['prior'] = ginfo.loc['filtered_n']/motif_all
    ginfo.loc['prior'] = ginfo.loc['motif_n'] / motif_all

    # normalized probability is raw count + prior / row_sum
    prob = motif.add(ginfo.loc['prior'])
    prob = prob.div(prob.sum(axis=1), axis=0)

    # Calculate pairwise Jensen - Shannon Distance
    # ensenshannon() calculates the square root of JS Divergence
    features = prob.index
    feature_n = len(features)
    js_distances = {}

    for f1, f2 in combinations(features, 2):
        p = prob.loc[f1]
        q = prob.loc[f2]
        # compute distance, then square to get divergence
        distance = jensenshannon(p, q)
        js_distances[(f1, f2)] = distance * distance  # Jensen-Shannon Divergence

    n = 0
    features_found = defaultdict(int)
    for pair, jsd in sorted(js_distances.items(), key=lambda p: p[1], reverse=True):
        n += 1
        pmax = motif.loc[pair[0]].idxmax()
        qmax = motif.loc[pair[1]].idxmax()
        features_found[pair[0]] += 1
        features_found[pair[1]] += 1
        feature_ratio = len(features_found) / feature_n
        print(f"{n}\t{pair[0]}\t{pair[1]}\t{jsd:.3f}\t{pmax}|{qmax}\t{len(features_found)}\t{feature_ratio:.3f}")
        if feature_ratio > target_ratio:
            break

    print(f'\nselected features({len(features_found)})')
    for motifid in sorted(features_found):
        # print(motif)
        # mprob = prob.loc[motif]
        n = motif.loc[motifid].sum()
        probstr = prob.loc[motifid].to_string(index=False, header=False, float_format='%7.3f')
        probstr = probstr.replace('\n', '\t')
        print(f'{motifid}\t{n}\t{probstr}')

    pd.options.display.max_rows = 2000
    pd.options.display.max_columns = 20
    pd.options.display.max_colwidth = 100
    pd.options.display.width = 500
    pd.set_option("display.precision", 3)

    print(f'\n{ginfo}')

    exit(0)
