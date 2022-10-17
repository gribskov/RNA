import glob
import os

fptdict = {}
# outerloop over fingerprint files

fptfiles = glob.glob('*.fpt')
for fpt in fptfiles:
    # fasta will be the filename with complete directory path

    # if you need just the filename without directory
    base = os.path.basename(fpt)
    dist = open(fpt, 'r')
    motifs = {}
    motifsfound = False
    # loop over motifs in fingerprint
    for line in dist:
        if line.find('- motif') > -1:
            motifsfound = True
            continue

        if motifsfound:
            field = line.split(':', maxsplit=1)
        else:
            continue

        motifs[field[0]] = int(field[1].rstrip())
        # print(motifs)

    fptdict[base] = motifs

motifcount = {}
for fpt in fptdict:
    for motif in fptdict[fpt]:
        if motif in motifcount:
            motifcount[motif]+=1
        else:
            motifcount[motif]=1
# print(motifcount)

for item in sorted(motifcount,key=lambda x:motifcount[x],reverse=True):
        print(item,motifcount[item])