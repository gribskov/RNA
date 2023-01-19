"""=================================================================================================
Run multiple jobs using subprocess.Popen().
Usage
    rna_manager --jobs <n> <base_directory>
    rna_manager --j <n> <base_directory>

--jobs is the number of processes to run simultaneously
see Pipeline:arg_get or run with -h to see additional arguments

Ben Iovino    19 April 2022
Michael Gribskov    3 May 2022

from github.com/biovino/biol494/rna_manager.py commit 7615726bdcc439741db9ebca8805da54b2b386c2
================================================================================================="""

import subprocess as sub
import time
import os
import glob
import sys
import argparse


class Pipeline():
    '''=============================================================================================
    =============================================================================================='''

    def __init__(self, base='./'):
        '''-----------------------------------------------------------------------------------------
        manager object holds locations of files and directories
        -----------------------------------------------------------------------------------------'''
        # define base here so it can be used for defaults in get_args()
        self.base = base
        if not self.base.endswith('/'):
            self.base += '/'

        # these defaults are overridden by command line arguments, if present. see arg_get()
        #
        self.log = './log/'
        self.fasta = './fasta/*.fa'
        self.ct = './ctfiles/'
        self.xios = './xiosfiles/'
        self.fpt = './fpt/'
        self.jobs = 20

        # instance variables
        self.managerlog = None
        self.errorlog = None

        self.stage = []
        self.source = ''  # directory for input to the current stage
        self.current = ''
        self.complete = {}

        self.jobid = 0
        self.running = 0
        self.total = 0
        self.started = 0
        self.finished = 0
        self.succeeded = 0
        self.failed = 0
        self.delay = 5
        self.joblist = []

        self.args = self.arg_get()
        args = self.args

    @staticmethod
    def arg_formatter(prog):
        """---------------------------------------------------------------------------------------------
        Set up formatting for help. Used in arg_get.

        :param prog:
        :return: argparse formatter class
        ---------------------------------------------------------------------------------------------"""
        return argparse.HelpFormatter(prog, max_help_position=60, width=120)

    def arg_get(self):
        """-----------------------------------------------------------------------------------------
        Command line arguments

        :return: command line namespace
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
        # paths to executables

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
        commandline.add_argument('-i', '--fasta', '--input',
                                 metavar='DIR_NAME',
                                 help='path to input fasta files FastA files (%(default)s)',
                                 default=f'{base}fasta/*.fa')
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

    def stageadd(self, description):
        """-----------------------------------------------------------------------------------------
        Append a stage description to self.stage
        Description is a dictionary with keys:
            stage - string, name of the stage
            command - string, the base command being run
            options - list of strings, command options
            source - list or strings, directory and suffix of input files
            dirs - list of strings, output directories (existence will be verified)

            example:
            {'stage':  'xios',
                       'command': f'python3 {workflow.python}/xios_from_rnastructure.py',
                       'options': ['-q ',
                                   f'-c {workflow.ct}',
                                   f'-x {workflow.xios}',
                                   f'-f $FILE',
                                   f'-y {workflow.python}',
                                   f'-r {workflow.RNAstructure}',
                                   ],
                       'source':  [workflow.fasta,'.fa'],
                       'dirs':    [workflow.ct, workflow.xios]
                       })


        :param description: dict
        :return: int, total number of stages
        -----------------------------------------------------------------------------------------"""
        self.stage.append(description)
        return len(self.stage)

    def fastforward(self):
        """-----------------------------------------------------------------------------------------
        Read the most recent logfile and make a list of completed files for each stage 
        (self.complete).
        Reinitializes self.complete

        :return: int - complete_n
        -----------------------------------------------------------------------------------------"""
        # initialize lists of completed files for each stage
        self.complete = {self.stage[i]['stage']: [] for i in range(len(self.stage))}

        # check the log directory and create if necessary
        create_files_and_stop = Pipeline.dircheck(self.log)

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
            # recent logfile exists, make a list of completed files
            complete_n = 0
            for line in manager:
                # print('manager line', line)
                info = self.logread(line, 'complete')
                # print('info', info)
                if info:
                    complete_n += 1
                    # print(info['message'])
                    start = info['message'].find('file:')
                    file = info['message'][start + 5:]
                    self.complete[info['stage']].append(file)
                    self.logwrite('manager', 'complete', info['stage'], f'previous log;file:{file}')

            manager.close()

        return complete_n

    @staticmethod
    def logtime():
        """-----------------------------------------------------------------------------------------
        Create a time string for use in logs year month day hour min sec
        concatenated

        :return: str - YYYYMoDyHrMnSc
        -----------------------------------------------------------------------------------------"""
        return time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))

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
            dirlist.sort()
            # print(f'opening {dirlist[-1]}, alt {dirlist[0]}')
            return open(dirlist[-2], 'r')
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
        log_message = f'{Pipeline.logtime()}\t{stage}\t{tag:12s}\t{message}\n'
        fp.write(log_message)
        fp.flush()

        return log_message

    def logread(self, line, filter):
        """-----------------------------------------------------------------------------------------
        read the indicated log until a matching tag is found. return a hash with the log information
        or if nothing is found, an empty dict (False)
        
        :param log: string - name of log (e.g., manager, error)
        :param filter: string: tag field must match
        :return: dict - keys: time, stage, tag, message
        -----------------------------------------------------------------------------------------"""
        info = {}
        field = line.rstrip().split('\t')
        # print(f'logread {field}', end='\t')
        # print(field)
        if field[2].startswith(filter):
            # print('match', end='')
            info = {'time': field[0], 'stage': field[1], 'tag': field[2], 'message': field[3]}
        # print(f'info {info}')

        return info

    @staticmethod
    def get_file_list(fileglob, filter='.fa'):
        """-----------------------------------------------------------------------------------------
        uses a wildcard filename to produe a list of files

        :param directory: string - full path of desired files i.e. /scratch/scholar/user/data/
        :param filter: string - file must contain this string
        :return filelist: list - file names i.e. [fasta1.fa, fasta2.fa, fasta3.fa]
        -----------------------------------------------------------------------------------------"""
        # Initialize a list of file names
        if not fileglob.endswith('*'):
            fileglob += '*'
        target = fileglob + filter
        filelist = glob.glob(target)

        if filter:
            for file in filelist:
                # shouldn't really need to filter, but kept it anyway
                if filter not in file:
                    # append file name to list if it matches filter
                    filelist.remove(file)

        return filelist

    @staticmethod
    def dircheck(path):
        """-----------------------------------------------------------------------------------------
        Create directory if it doesn't exist already

        :param dirname: string - directory path
        :return: boolean - True if new directory created False otherwise
        -----------------------------------------------------------------------------------------"""
        if not os.path.isdir(path):
            os.mkdir(path)
            return True

        return False

    def manager(self):
        """-----------------------------------------------------------------------------------------
        Runs each stage keeping self.jobs running at a time.  Delay between polling is set by 
        self.delay
        -----------------------------------------------------------------------------------------"""
        # fastforward checks the log for files that have already been completed
        self.fastforward()
        files = []
        for stage in self.stage:
            s = stage['stage']
            self.source = stage['source']
            files.append(f'{s}:{len(self.complete[s])}')
        self.logwrite('manager', 'fastforward', f'init', ';'.join(files))

        for stage in self.stage:
            self.current = stage['stage']
            sourcedir, sourcetype = stage['source']
            # stagedir = getattr(self, self.source)
            # filelist = Pipeline.get_file_list(self.current)
            print(f'\nstage={stage["stage"]}\tsource={sourcedir}*{sourcetype}')
            filelist = Pipeline.get_file_list(sourcedir, sourcetype)

            # total number of jobs for this stage
            total = len(filelist)
            self.logwrite('manager', 'files', f'{stage["stage"]}', f'{total} files in {self.source}')

            # check that output directories exist
            for dir in stage['dirs']:
                Pipeline.dircheck(dir)

            while self.manager_startjobs(filelist, stage) or self.joblist:
                self.manager_polljobs(stage)

        return self.finished

    def manager_startjobs(self, filelist, stage):
        """-----------------------------------------------------------------------------------------
        process all files in filelist using the current stage command. run up to self.jobs at a time

        usage:
            while self.manager_startjobs(filelist, stage):
                self.manager_polljobs(stage)

        :param filelist: list - list of input files to process
        :param stage: dict - information describing current pipeline stage
        :return: boolean - True indicates there are more files to process
        -----------------------------------------------------------------------------------------"""
        this_stage = stage['stage']
        while filelist:
            # print(f'startjobs filelist:{len(filelist)}')

            if self.running < self.jobs:
                file = filelist.pop()
                if file in self.complete[this_stage]:
                    print(f'skipping {file}')
                    continue
                # print(f'loop filelist:{len(filelist)}')

                self.jobid += 1
                print(f'starting job{self.jobid} {file}')

                thiscommand = [stage['command']] + stage['options']
                thiscommand = [clause.replace('$FILE', file) for clause in thiscommand]

                self.logwrite('manager', 'start', stage['stage'], f'jobid:{self.jobid}; input:{file} ')
                self.logwrite('manager', 'command', stage['stage'], ' '.join(thiscommand))
                job = sub.Popen(' '.join(thiscommand), shell=True, stdout=self.errorlog, stderr=self.errorlog)
                self.joblist.append([self.jobid, job, file])
                self.running += 1

            else:
                # desired number of jobs are running so continue to polling
                break

        if filelist:
            # files remain to be processed
            return True
        else:
            return False

    def manager_polljobs(self, stage):
        """-----------------------------------------------------------------------------------------
        Poll the currently running jobs, and remove completed jobs from the joblist

        :param stage: dict - information describing current pipeline stage
        :return:
        -----------------------------------------------------------------------------------------"""
        # poll all jobs in joblist
        time.sleep(self.delay)
        to_remove = []
        for j in self.joblist:
            id, job, file = j
            # print(f'\tjob {id} ...', end='')
            result = job.poll()
            if result == None:
                # None indicates job is still running
                # print(f'job {id} still running')
                pass

            else:
                # job finished
                # print(f'job{id} finished')
                self.running -= 1
                self.finished += 1
                if result == 0:
                    # success
                    print(f'job{id} succeeded {file}')
                    self.logwrite('manager', 'complete', stage['stage'], f'jobid:{id};file:{file}')
                    self.succeeded += 1

                else:
                    # error
                    print(f'job{id} failed {file}')
                    self.logwrite('manager', 'fail', self.current, f'status:{result};file:{file}')
                    self.failed += 1

                # include the result in the remove list, it can't be removed here because it
                # changes the list (self.joblist) that is iterating
                to_remove.append(j)

        # remove all finished jobs. Couldn't do it above because it shortens the joblist and 
        # some jobs don't get polled
        for j in to_remove:
            self.joblist.remove(j)

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

workflow = Pipeline()
workflow.source = 'fasta'
workflow.delay = 5

# add commands
# TODO change this to be read from external file (probably YAML)
workflow.stage.append({'stage': 'xios',
                       'command': f'python3 {workflow.python}/xios_from_rnastructure.py',
                       'options': ['-q ',
                                   f'-c {workflow.ct}',
                                   f'-x {workflow.xios}',
                                   f'-f $FILE',
                                   f'-y {workflow.python}',
                                   f'-r {workflow.RNAstructure}',
                                   f'-d4,6 -w9,11',
                                   ],
                       'source': [workflow.fasta, '.fa'],
                       'dirs': [workflow.ct, workflow.xios]
                       })
workflow.stage.append({'stage': 'fpt',
                       'command': f'python3 {workflow.python}/fingerprint_random.py',
                       'options': [f'-r $FILE',
                                   f'-f auto',
                                   f'--motifdb {workflow.python + "/data/2to7stem.mdb.pkl"}',
                                   f'--outputdir {workflow.fpt}',
                                   f'-c 5 -l 100000',
                                   f'-q',
                                   f'-n'],
                       'source': [workflow.xios, '.xios'],
                       'dirs': [workflow.fpt]
                       })

workflow.manager()
workflow.cleanup()

exit(0)
