dist=open(tmRNA.Magnetococcus_marinus.CP000471.fpt,'r')
for line in 'dist':
    if line.find('-motif'):
        motifsfound=True
    if motifsfound:
        field=line.split(':',maxsplit=1)
        id=field[0].replace('"',")


