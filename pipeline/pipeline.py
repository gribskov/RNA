"""=================================================================================================
Run multiple jobs using subprocess.Popen()
Poll jobs until all jobs are done
================================================================================================="""
import subprocess as sub
import random
from time import sleep

jobs = []
n = 3

# logfile for output
log = open('sleep.log', 'wb')

# start n jobs, each job sleeps a random number of seconds, then terminates
for i in range(n):
    sec = random.randrange(15)
    command = ['py', 'sleep.py', f'{sec}']
    print(f'job {i}: {sec} seconds')
    # job = sub.Popen(command, shell=True, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
    # job = sub.Popen(command, shell=True, stdout=log, stderr=log)
    job = sub.Popen(command, shell=True)
    jobs.append(job)

# poll until all jobs finish
done = 0
delay = 2  # number of seconds to wait between polls
while done < n:
    print('\nPolling')
    for i in range(n):
        if jobs[i] == 'Done':
            continue

        print('    job {} ...'.format(i), end='')
        result = jobs[i].poll()
        # out, err = jobs[i].communicate()


        if result != None:
            print('finished')
            jobs[i] = 'Done'
            done += 1

        else:
            print('still running')
            # print('job {}:result={}'.format(i, result))

    sleep(delay)

log.flush()
log.close()
