from topology import Topology, PairRNA


if __name__ == '__main__':

    topo = PairRNA(inlist=[0,0])
    print(topo)

    level_max = 3
    level = 0
    while level < level_max:
        level += 1
        for left in range(len(topo)):
            for right in range(left+1, len(topo)+1):
                topo.pairs = [[left,right]] + topo.pairs
                for i in range(1,len(topo.pairs)):
                    for j in (1,2):
                        if topo.pairs[i][j] >= left:
                            topo.pairs[i][j] += 1
                        if topo.pairs[i][j] >= right:
                            topo.pairs[i][j] += 2
                print(topo)





    exit(0)
