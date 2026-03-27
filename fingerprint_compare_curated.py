#===================================================================================================
# compare fingerprints to existing curated results
#
# mgribsko      3/10/26
#===================================================================================================
import sys
import os
import glob
from lxml import etree
from fingerprint import FingerprintSet, Fingerprint

if __name__ == '__main__':
    curdir = f'{sys.argv[1]}/*.fpt'
    testparent = f'{sys.argv[2]}'


    # read curated fingerprints
    fpt_list = glob.glob(curdir)
    curated = {}
    count = 0
    for this_fpt in fpt_list:
        f = Fingerprint()
        f.readYAML(this_fpt)
        curated[f.information['File']] = f
        count += 1
        print(f'\t{count}\t{f.information['File']}')

    print(f'{len(curated)} curated fingerprints read from {curdir}')
    
    # read and compare test fingerprints
    # fpt_list = glob.glob(testdir)
    # test = {}
    # count = 0
    # for this_fpt in fpt_list:
    #     f = Fingerprint()
    #     f.readYAML(this_fpt)
    #     name = os.path.basename(f.information['File'])
    #     test[name] = f
    #     count += 1
    #     print(f'\t{count}\t{name}')

    #print(f'{len(test)} test fingerprints read from {testparent}')
    # all directories in testparent
    out =open('summary.out', 'w')
    directories = [entry.name for entry in os.scandir(testparent) if entry.is_dir()]
    dirlist = []
    fptlist = {}
    dcount = 0
    fptcount = 0
    for dir in directories:
        fpt = f'{testparent}/{dir}/fpt'
        #print(f'{dcount}: {fpt}')
        dcount += 1
        with os.scandir(fpt) as entries:
            for entry in entries:
                #print(f'{entry.name}')
                if entry.name.endswith('.fpt'):
                    # each file in the fpt directory
                    fptcount += 1
                    t = Fingerprint()
                    t.readYAML(f'{fpt}/{entry.name}')
                    try:
                        c = curated[entry.name]
                        #t = test[fpt]
                    except KeyError as e:
                        print(f'!Error: {fpt} not found in test')
                        continue
                    clist = []
                    tlist = []
                    both = []
                    for m in c.motif:
                        if m in t.motif:
                            both.append(m)
                        else:
                            clist.append(m)
                    for m in t.motif:
                        if m not in c.motif:
                            tlist.append(m)

                    nb = len(both)
                    nc = len(clist)
                    nt = len(tlist)
                    recall = nb / (nb + nc)
                    precision = nb / (nb + nt)
                    try:
                        jaccard = nb / (nt + nc + nb)
                    except ZeroDivisionError:
                        jaccard = 0
                    out.write(f'{recall:6.3f}  {precision:6.3f} {jaccard:6.3f}  ')
                    out.write(f'{len(clist)}\t{len(both)}\t{len(tlist)}\t{dir}\t{entry.name}\n')


            fptlist[dir] = fptcount
    # for fpt in sorted(curated):
    #     try:
    #         c = curated[fpt]
    #         t = test[fpt]
    #     except KeyError as e:
    #         print(f'{fpt} not found in test')
    #         continue
    #     clist = []
    #     tlist = []
    #     both = []
    #     for  m in c.motif:
    #         if m in t.motif:
    #             both.append(m)
    #         else:
    #             clist.append(m)
    #     for m in t.motif:
    #         if m not in c.motif:
    #             tlist.append(m)
    #
    #     nb = len(both)
    #     nc = len(clist)
    #     nt = len(tlist)
    #     recall = nb / (nb + nc)
    #     precision = nb / (nb + nt)
    #     try:
    #         jaccard = nb / (nt + nc + nb)
    #     except ZeroDivisionError:
    #         jaccard = 0
    #     print(f'{recall:6.3f}  {precision:6.3f} {jaccard:6.3f}  ', end='')
    #     print(f'{len(clist)}\t{len(both)}\t{len(tlist)}\t{fpt}')


    exit(0)
