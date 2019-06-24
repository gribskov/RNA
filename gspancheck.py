"""=================================================================================================
more unique graphs are found when forward and backward structures are porocessed separately,
implying that gspan givws different answers.  this programs genreates forward and backaward
versions of the 7 stem motif list

Michael Gribskov 21 June 2019
================================================================================================="""
from motifdb import SerialRNA
from xios import Xios
from gspan import Gspan

field = []
n = 0
motif = open('data/stem12.list')

for line in motif:
    n += 1
    line = motif.readline()
    field = line.split()
    # print('{} {}'.format(n, field[1]))

    forward = SerialRNA(field[1])
    forward.canonical()
    backward = SerialRNA(forward)
    backward.reverse()
    backward.canonical()

    # forward min dfs
    fxios = Xios(serial=forward)
    fgspan = Gspan(fxios)
    fdfs = fgspan.minDFS()
    fdfshex = fdfs.ascii_encode()

    # backward min dfs
    bxios = Xios(serial=backward)
    bgspan = Gspan(bxios)
    bdfs = bgspan.minDFS()
    bdfshex = bdfs.ascii_encode()

    if bdfshex != fdfshex:
        print('{}\t {}\t{}'.format(n, fxios, bxios))
        print('\t\tforward {}    backward{}'.format(forward, backward))
        print('\t\t{}\t{}\n\t\t{}\t{}'.format(fdfs, fgspan.g2d, bdfs, bgspan.g2d))