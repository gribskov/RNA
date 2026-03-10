import os

# cancel slurm jobs in todelete.list
# generate the list with squeue: squeue -A gcore > todelete.list 
#
td = open('todelete.list', 'r' )
for line in td:
    field = line.split()
    print(field[0])
    cancel =f'scancel {field[0]}'
    os.system(cancel)
