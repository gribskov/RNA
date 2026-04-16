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



    exit(0)
