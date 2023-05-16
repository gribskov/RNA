"""
test script for windows, just sleeps for a random number of seconds between 0 and 30
"""
import random
import sys
import time

rtime= random.random() * 60
start = time.time()
sys.stdout.write(f'out:sleeping {rtime:.1f} seconds, job:{sys.argv[1]}\tstarted:{start:.1f}\n')
time.sleep(rtime)
sys.stderr.write(f'\terr:sleeping {rtime:.1f} seconds, job:{sys.argv[1]}\tstarted:{start:.1f}\n')

exit(0)