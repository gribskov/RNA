"""=================================================================================================
reduce dimensionality of motif features?
================================================================================================="""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn import datasets

from fingerprint import FingerprintMatrix
import os

if __name__ == '__main__':

    np.random.seed(9237)
    # iris = datasets.load_iris()
    # X = iris.data
    # y = iris.target
    #
    # fig = plt.figure(1, figsize=(4, 3))
    # plt.clf()

    # ax = fig.add_subplot(111, projection="3d", elev=48, azim=134)
    # ax.set_position([0, 0, 0.95, 1])
    #
    # plt.cla()
    # pca = decomposition.PCA(n_components=3)
    # pca.fit(X)
    # X = pca.transform(X)
    #
    # for name, label in [("Setosa", 0), ("Versicolour", 1), ("Virginica", 2)]:
    #     ax.text3D(
    #         X[y == label, 0].mean(),
    #         X[y == label, 1].mean() + 1.5,
    #         X[y == label, 2].mean(),
    #         name,
    #         horizontalalignment="center",
    #         bbox=dict(alpha=0.5, edgecolor="w", facecolor="w"),
    #     )
    # # Reorder the labels to have colors matching the cluster results
    # y = np.choose(y, [1, 2, 0]).astype(float)
    # ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap=plt.cm.nipy_spectral, edgecolor="k")
    #
    # ax.xaxis.set_ticklabels([])
    # ax.yaxis.set_ticklabels([])
    # ax.zaxis.set_ticklabels([])
    #
    # plt.show()

    # file = 'fmatrix.tsv'
    # read from pickled fingerprint matrix with fingerprints = rows and motifs = columns
    picklefile = 'fmatrix.pkl'
    fmat = FingerprintMatrix.unpickle(picklefile)
    print(f'{len(fmat.fpt)} fingerprints with {len(fmat.motifs)} motifs unpickled')
    motif_mincount = 10
    motif_maxcount = 140
    fmat.select_min_max(motif_mincount, motif_maxcount, False, recalculate=True)
    f = fmat.fpt
    motif = np.array([[f[row][col] for col in range(len(f[0]))] for row in range(len(f))])
    colnames = fmat.motifs_selected()
    print(f'{len(colnames)} motifs selected')

    # get fingerprint names and trim
    rownamedict = fmat.fpt_id
    rownames = []
    for name in rownamedict:
        name = name.replace('.fpt','')
        name = name.replace('.xios.out', '')
        name = os.path.basename(name)
        rownames.append(name)

    # groups/labels
    # for the curated data, if the first token after splitting on . and _ matches, it is the same group
    groups = {}
    gindex = []
    g = []
    for name in rownames:
        token = name.split('.')
        token = token[0].split('_')[0]
        if token not in groups:
            groups[token] = len(groups)
            gindex.append(token)

        g.append(groups[token])
    print(f'{len(groups)} groups identified')
    component_n = 3
    print(f'starting PCA (components={component_n})...')
    pca = decomposition.PCA(n_components=component_n)
    pca.fit(motif)
    X = pca.transform(motif)
    components = pca.components_
    print(f'PCA complete')
    for c in pca.components_:
        count = 0
        for i in sorted(range(len(c)), key=lambda i: c[i], reverse=True):
            print(f'{i}\t{c[i]:.3g}\t{rownames[i]}')
            count += 1
            if count > 5:
                break
        print()
    print(pca.explained_variance_ratio_)
    print(pca.singular_values_)
    fig = plt.figure(1, figsize=(4, 3))
    plt.clf()

    ax = fig.add_subplot(111, projection="3d", elev=48, azim=134)
    ax.set_position([0, 0, 0.95, 1])
    # y = np.choose(y, [1, 2, 0]).astype(float)
    # scatter = ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=g, cmap=plt.cm.nipy_spectral, label=g, edgecolor="k")
    scatter = ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=g, cmap=plt.cm.Accent, label=g, edgecolor="k", s=50, alpha = 0.8)
    a = scatter.legend_elements()

    for i in range(len(a[1])):
        a[1][i] = gindex[i]

    legend1 = ax.legend(*a,
                        loc="upper left", title="Groups")
    # gmean = [0 for _ in g]
    # for name, label in groups.items():
    #    gmean[]
    # ax.text3D(X[g==label, 0].mean(),
    #     name,
    #     horizontalalignment="center",
    #     bbox=dict(alpha=0.5, edgecolor="w", facecolor="w"),
    # )
    plt.show()
    exit(0)
