"""=================================================================================================
Run multiple jobs using subprocess.Popen().
total is the total number of jobs to run.
run is the number to keep running at the same time
Ben Iovino    19 April 2022

from github.com/biovino/biol494/rna_manager.py commit 7615726bdcc439741db9ebca8805da54b2b386c2
================================================================================================="""

import subprocess as sub
import time
import os
import sys
import argparse


class Pipeline():
    '''=============================================================================================
    =============================================================================================='''

    def __init__(self, base='./'):
        '''-----------------------------------------------------------------------------------------
        manager object holds locations of files and directories
        -----------------------------------------------------------------------------------------'''
        self.jobs = 20
        self.running = 0
        self.total = 0
        self.finished = 0
        self.delay = 5
        self.joblist = []

        self.base = base
        if not self.base.endswith('/'):
            self.base += '/'

        self.RNAastructure = '../RNAstructure'
        self.python = '../RNA'

        self.log = self.base + 'log/'
        self.fasta = self.base + 'fasta/'
        self.ct = self.base + 'ct/'
        self.xios = self.base + 'xios/'
        self.fpt = self.base + 'fpt/'

        self.args = self.arg_get()

        self.managerlog = None
        self.errorlog = None

        self.stage = []
        self.source = ''  # directory for input to the current stage
        self.current = ''
        self.complete = {}

    @staticmethod
    def arg_formatter(prog):
        """---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Set up formatting for help
        :param prog:
        :return: argparse formatter class
        ---------------------------------------------------------------------------------------------"""
        return argparse.HelpFormatter(prog, max_help_position=60, width=120)

    def arg_get(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        base = self.base
        if not base.endswith('/'):
            base += '/'

        commandline = argparse.ArgumentParser(
            description='Run XIOS fingerprint pipeline', formatter_class=Pipeline.arg_formatter)

        commandline.add_argument('-j', '--jobs',
                                 help='number of concurrent jobs to run',
                                 type=int,
                                 default=20)
        commandline.add_argument('-q', '--quiet',
                                 help='run with minimal output to terminal',
                                 action='store_true')
        commandline.add_argument('-p', '--passthrough',
                                 metavar='COMMAND_STRING',
                                 help='command options to pass through to pipeline executables (%(default)s)',
                                 default='NONE')
        # # paths to executables

        commandline.add_argument('-R', '--RNAstructure',
                                 metavar='DIRECTORY',
                                 help='directory with RNAstructure package executables (%(default)s)',
                                 default='../RNAstructure')
        commandline.add_argument('-P', '--python',
                                 metavar='DIRECTORY',
                                 help='directory with XIOS python executable3s (%(default)s)',
                                 default='../RNA')

        # alternate names for input and output directories
        commandline.add_argument('-l', '--log',
                                 metavar='DIR_NAME',
                                 help='log directory (%(default)s)',
                                 default=f'{base}log/')
        commandline.add_argument('-i', '--fasta',
                                 metavar='DIR_NAME',
                                 help='subdirectory with initial FastA files (%(default)s)',
                                 default=f'{base}fasta/')
        commandline.add_argument('-c', '--ctdir',
                                 metavar='DIR_NAME',
                                 help='subdirectory for CT files (%(default)s)',
                                 default=f'{base}ct/')
        commandline.add_argument('-x', '--xiosdir',
                                 metavar='DIR_NAME',
                                 help='subdirectory for XIOS files (%(default)s)',
                                 default=f'{base}xios/')
        commandline.add_argument('-f', '--fingerprint',
                                 metavar='DIR_NAME',
                                 help='subdirectory for fingerprint files (%(default)s)',
                                 default=f'{base}fpt/')

        commandline.add_argument('base', type=str,
                                 help='Base directory for project (%(default)s)',
                                 default=f'{base}')

        args = commandline.parse_args()
        self.RNAstructure = args.RNAstructure
        self.python = args.python

        # change object values if specified on command line
        if not args.base.endswith('/'):
            args.base += '/'
        self.base = args.base if args.base else self.base
        self.log = args.log if args.log else self.log
        self.fasta = args.fasta if args.fasta else self.fasta
        self.ct = args.ctdir if args.ctdir else self.ct
        self.xios = args.xiosdir if args.xiosdir else self.xios
        self.fpt = args.fingerprint if args.fingerprint else self.fpt

        # commandline.add_argument('-p', '--percent',
        #                          help='Fold: percent cutoff for suboptimal structures (%(default)s)',
        #                          default=10)
        # commandline.add_argument('-m', '--maximum',
        #                          help='Fold: maximum number of suboptimal structures (%(default)s)',
        #                          default=100)
        # commandline.add_argument('-w', '--window',
        #                          help='Fold: window for suboptimal structures (%(default)s) or comma separated',
        #                          default='4')
        # commandline.add_argument('-d', '--ddG',
        #                          help='Mergestems: delta deltaG limit for supoptimal structures (%(default)s) '
        #                               'or comma separated',
        #                          default='5')

        return args

    def stageadd(self, name, description):
        """-----------------------------------------------------------------------------------------
        Append a stage description to self.stage
        TODO add details of description

        :param description: dict
        :return: int, total number of stages
        -----------------------------------------------------------------------------------------"""
        self.stage.append(description)
        return len(self.stage)

    # def stagefind(self):
    #     """-----------------------------------------------------------------------------------------
    #     Find the current stage by examining whether the stage directory has been created.  The
    #     output of each stage goes in a directory with the same name
    #
    #     :return: string, stage name
    #     -----------------------------------------------------------------------------------------"""
    #     for stage in self.stage:
    #         path = self.base + stage['stage']
    #         if os.path.isdir(path):
    #             # current will be the last one matched
    #             self.source = stage
    #             self.current = stage
    #
    #     if not self.current:
    #         # first stage if none have been run before
    #         self.current = self.stage[0]['stage']
    #
    #     return self.current

    def fastforward(self):
        """-----------------------------------------------------------------------------------------
        Read the most recent logfile and make a list of completed files for each stage.
        Reinitializes self.complete

        :return: complete_n
        -----------------------------------------------------------------------------------------"""
        # initialize lists of completed files for each stage
        self.complete = {self.stage[i]['stage']: [] for i in range(len(self.stage))}

        # check the log directory and create if necessary
        create_files_and_stop = False
        if not os.path.isdir(self.log):
            # if directory doesn't exist, create it.
            os.mkdir(self.log)
            create_files_and_stop = True

        # create log files for the current run
        timestamp = Pipeline.logtime()
        self.managerlog = open(f'{self.log}{timestamp}_manager.log', 'w')
        self.errorlog = open(f'{self.log}{timestamp}_error.log', 'w')
        self.logwrite('manager', 'start_run', 'init', 'initialize logs')

        if create_files_and_stop:
            return 0

        # if there are existing logs, identify the most recent one and read the
        # completed files
        manager = self.logrecent('manager')
        if manager:
            # recent logfile exists, read in completed files
            complete_n = 0
            for line in manager:
                info = self.logread(manager, 'complete')
                if info:
                    complete_n += 1
                    self.complete[info['stage']].append(info['message'])
                    self.logwrite('manager', 'complete', info['stage'], 'read from previous log')

            manager.close()

        return complete_n

    @staticmethod
    def logtime():
        """
        Create a time string for use in logs year month day hour min sec
        concatenated

        :return:
        """
        t = time.localtime()
        return time.strftime('%Y%m%d%H%M%S')

    def logrecent(self, logname):
        """-----------------------------------------------------------------------------------------
        Find and open the most recent manager log for reading

        :return: file pointer of open file
        -----------------------------------------------------------------------------------------"""
        dirlist = []
        for file in os.listdir(self.log):
            d = os.path.join(self.log, file)
            if logname in d:
                dirlist.append(d)

        if dirlist:
            return open(dirlist[-1], 'r')
        else:
            return None

    def logwrite(self, log, tag, stage, message):
        """-----------------------------------------------------------------------------------------
        add a timestamped message to the named log

        :param log: string - name of log (e.g., manager, error)
        :param tag: string - tag describing action
        :param stage: string - pipeline stage
        :param message: string - message
        :return: string - message written to log
        -----------------------------------------------------------------------------------------"""
        fp = None
        if log == 'manager':
            fp = self.managerlog
        elif log == 'error':
            fp = self.errorlog
        else:
            self.logwrite('error', 'unknown_logfile', stage, '')

        # log_message = f'{tag}\t{stage}\t{message}\t{self.logtime()}\n'
        log_message = f'{self.logtime()}\t{stage}\t{tag:12s}\t{message}\n'
        fp.write(log_message)
        fp.flush()

        return log_message

    def logread(self, fp, filter):
        """-----------------------------------------------------------------------------------------
        read the indicated log until a matching tag is found. return a hash with the log information
        or if nothing is found, an empty dict (False)
        
        :param log: string - name of log (e.g., manager, error)
        :param filter: string: tag field must match
        :return: dict - keys: time, stage, tag, message
        -----------------------------------------------------------------------------------------"""
        info = {}
        for line in fp:
            field = line.rstrip().split('\t')
            if info['tag'] == 'completed':
                info = {'time':field[0], 'stage':field[1], 'tag':field[2], 'message':field[3]}
                break

        return info

    @staticmethod
    def get_file_list(directory, filter='.fa'):
        """-----------------------------------------------------------------------------------------
        Goes through the input directory and returns a list of each file in the directory. Only
        files that contain the string filter are selected

        :param directory: string - full path of desired files i.e. /scratch/scholar/user/data/
        :param filter: string - file must contain this string
        :return filelist: list - file names i.e. [fasta1.fa, fasta2.fa, fasta3.fa]
        -----------------------------------------------------------------------------------------"""
        # Initialize a list of file names
        filelist = list()

        for file in os.listdir(directory):
            # Read each file name in the directory path
            if filter in file:
                # append file name to list if it matches filter
                filelist.append(directory+file)

        # filelist.sort()   # probably better to not sort

        return filelist

    def manager(self):
        """-----------------------------------------------------------------------------------------
        Accepts a directory and list of files, runs xios_from_rnastructure.py on each file. Runs a
        certain amount of jobs at once, logs xios_from_rnastructure.py output and logs fasta files
        that went through complete job. Polls every x amount of seconds to check which jobs are
        complete.
        :param base: directory of directories i.e. /scratch/scholar/user/data/
        :param pythonexe: directory with python executable files
        :param rnaexe: directory with RNAstructure executable files
        :param jobs: number of jobs to run at once
        :param w: window for RNAstructure
        :param d: delta delta G for RNAstructure
        :param filelist: list of files i.e. [fasta1.fa, fasta2.fa, fasta3.fa]
        :param startwith: index of filelist to start running
        :param directory: directory of fasta files i.e. /scratch/scholar/user/data/avocado
        -----------------------------------------------------------------------------------------"""
        # self.logwrite('manager', 'start', 'init', '')

        self.fastforward()
        files = []
        for stage in self.stage:
            s = stage['stage']
            files.append(f'{s}:{len(self.complete[s])}')
        self.logwrite('manager', 'fastforward', f'init', ','.join(files))

        for stage in self.stage:
            stagedir = getattr(self, self.source)
            filelist = Pipeline.get_file_list(stagedir)

            # total number of jobs for this stage
            total = len(filelist)
            self.logwrite('manager', 'files', f'{stage["stage"]}', f'{total} files in {self.source}')

            # TODO check that output directory exists

            while self.manager_startjobs(filelist, stage):
                self.manager_polljobs()

        return self.finished

    def manager_startjobs(self, filelist, stage):
        """-----------------------------------------------------------------------------------------
        start n jobs, each job sleeps a random number of seconds, then terminates

        :param filelist:
        :return:
        -----------------------------------------------------------------------------------------"""
        job_id = 0

        while filelist:
            file = filelist.pop()
            if file in self.complete:
                continue
            if self.running < self.jobs:
                job_id += 1
                # fasta = filelist[startwith + total_started]
                # command = f'python {pythonexe}/xios_from_rnastructure.py -i {directory} ' \
                #           f'-c {directory}/ctfiles -x {directory}/xiosfiles -f {fasta} -y {pythonexe} -r {rnaexe}'

                command = stage['command']
                self.logwrite('manager', 'start', stage['stage'], f'jobid:{job_id}; input:{file} ')
                # job = sub.Popen(command, shell=True, stdout=xios_log, stderr=xios_log)
                job = 'test'
                self.joblist.append(job_id)
                self.running += 1
                # total_started += 1
            else:
                break

        if filelist:
            # files remain to be processed
            return True
        else:
            return False

    def manager_polljobs(self):
        """-----------------------------------------------------------------------------------------
        Poll the currently running jobs, and remove completed jobs from the joblist

        :return:
        -----------------------------------------------------------------------------------------"""
        # poll all jobs in joblist
        time.sleep(self.delay)
        to_remove = []
        for j in self.joblist:
            id, job = j
            # print(f'\tjob {id} ...', end='')
            result = job.poll()
            if result != None:
                # None indicates job is still running
                j.append(result)
                print(f'finished')
                to_remove.append(j)

            else:
                print('still running')

        # remove all finished jobs. Can't do it above because it shortens the joblist and some then
        # some jobs don't get polled
        for j in to_remove:

            # check for exit status
            # write fasta file name to file, and other marked attributes
            if j[3] == 2:
                # error
                self.logwrite('error', 'fail', self.current, f'status:{j[3]}')
            else:
                # success
                self.logwrite('manager', 'completed', self.current, f'jobid:{j[0]}')

            # remove finished jobs
            self.joblist.remove(j)
            self.running -= 1
            self.finished += 1

    def cleanup(self):
        """-----------------------------------------------------------------------------------------
        End of run cleanup.
        Close files

        :return: True
        -----------------------------------------------------------------------------------------"""
        self.managerlog.close()
        self.errorlog.close()

        return True


# --------------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------------

"""-------------------------------------------------------------------------------------------------
base is the main project directory, each step of the pipeline creates a directory under base
    ct
    xios
    fpt
    
Initializes directory, gets list of files from get_file_list(), checks if there is a log file and
reads it to find last fasta file worked on, and sends to manager()
-------------------------------------------------------------------------------------------------"""
# TODO add to passthrough arguments? or create command line switches?
# w = int(sys.argv[5])  # window param for xios_from_rnastructure.py
# d = int(sys.argv[6])  # delta delta G param for xios_from_rnastructure.py

workflow = Pipeline(base='data/')
# workflow.check_directory('log')
#
# command = f'python {pythonexe}/xios_from_rnastructure.py -i {directory} ' \
#           f'-c {directory}/ctfiles -x {directory}/xiosfiles -f {fasta} -y {pythonexe} -r {rnaexe}'
# add commands
workflow.stage.append({'stage':   'xios',
                       'command': f'python {workflow.python}/xios_from_rnastructure.py',
                       'options': [f'-c {workflow.ct}',
                                   f'-x {workflow.xios}',
                                   f'-f {workflow.fasta}',
                                   f'-y {workflow.python}',
                                   f'-r (workflow.RNAstructure',
                                   ]
                       })
workflow.source = 'fasta'
workflow.manager()
workflow.cleanup()

exit(0)
