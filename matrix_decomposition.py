"""=================================================================================================
reduce dimensionality of motif features?
================================================================================================="""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn import datasets


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
    file = 'fmatrix.tsv'
    # motif = pd.read_csv(file, sep='\t', index_col=0)
    # motif = np.loadtxt(file, dtype=int, usecols=(1,2,3,4,5,6,7,8), skiprows=1)
    motif = np.loadtxt(file, dtype=str, usecols=(1,2,3,4,5,6,7,8))
    colnames = motif[0,:]
    motif = np.loadtxt(file, dtype=str, skiprows=1)
    rownames = motif[:,0]
    values = np.transpose(motif[:,1:9].astype(np.int32))

    # print(motif.head())
    # print(motif.columns)
    # print(motif.info())
    # groups/labels
    g = []
    for col in colnames:
        print(col)
        if col.find('5S')>-1:
            g.append(0)
        elif col.find('tmRNA')>-1:
            g.append(1)
        elif col.find('rnasep')>-1:
            g.append(2)



    pca = decomposition.PCA(n_components=3)
    pca.fit(values)
    X = pca.transform(values)
    components = pca.components_
    for c in pca.components_:
        count = 0
        for i in sorted(range(len(c)), key=lambda i:c[i], reverse=True):
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
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=g, cmap=plt.cm.nipy_spectral, edgecolor="k")
    plt.show()
    exit(0)


