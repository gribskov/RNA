"""=================================================================================================
A simple program that sleeps for a number of seconds
use for testing popen submission

usage
    sleep n

    n = number of seconds to sleep

returns start and stop times and a unique ID

Michael Gribskov     03 March 2022
================================================================================================="""
import sys
import time
import uuid

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    u = uuid.uuid1()
    uid = u.hex[:5]

    t = time.localtime()
    now = time.strftime("%I:%M:%S ", t)
    sys.stdout.write(f'{uid} started {now}\n')

    time.sleep(int(sys.argv[1]))

    t = time.localtime()
    now = time.strftime("%I:%M:%S ", t)
    sys.stdout.write(f'{uid} finished {now}\n')

    exit(0)
