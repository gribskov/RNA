import os.path
import subprocess as sub
import shutil
import sys
import time
import argparse
import yaml


def arg_formatter(prog):
    """---------------------------------------------------------------------------------------------
    Set up formatting for help. Used in arg_get.

    :param prog:
    :return: argparse formatter class
    ---------------------------------------------------------------------------------------------"""
    return argparse.HelpFormatter(prog, max_help_position=60, width=120)


def getoptions():
    """-----------------------------------------------------------------------------------------
    Command line arguments:
    project: project directory, all stages will be in subdirectories

    r, denovo      remove previous results and perform all stages
    c, continue     continue existing run, skipping commands that are marked as complete

    options:
    j, jobs         number of concurrent jobs to run
    l, log          directory for log files

    :return: command line namespace
    -----------------------------------------------------------------------------------------"""
    project_default = './'
    commandline = argparse.ArgumentParser(
        description='Run workflow', formatter_class=arg_formatter)

    commandline.add_argument('project', type=str,
                             help='project directory for project (%(default)s)',
                             default=project_default)

    commandline.add_argument('-r', '--denovo',
                             help='Erase directories/files and denovo (%(default)s)',
                             default=False, action='store_true')

    commandline.add_argument('-w', '--workflow', type=str,
                             help='YAML workflow plan for project (%(default)s)',
                             default=f'workflow.yaml')

    commandline.add_argument('-j', '--jobs',
                             help='number of concurrent jobs to run',
                             type=int,
                             default=20)

    commandline.add_argument('-q', '--quiet',
                             help='run with minimal output to terminal',
                             action='store_true')

    commandline.add_argument('-l', '--log',
                             metavar='DIR_NAME',
                             help='log directory (%(default)s)',
                             default=f'{project_default}/log/')

    args = commandline.parse_args()

    # project directory must not end in /
    if not args.project.endswith('/'):
        args.project = args.project.rstrip('/')

    return vars(args)  # convert namespace to dict


class Workflow:
    """=============================================================================================
    python manager.py

    The workflow manager is responsible for reading the workflow plan and converting it to a series
    of executable commands. The files of command can then be executed in parallel by Executor.

    Workflow runs in two modes, controlled by --denovo:
    denovo: delete any existing directories/files and start a new project from scratch
    fastforward: (default) use existing directories and files. Execute only stages and commands
        that haven't completed

    main setup:
        create project directory, project is the base name for the run and for the main directory,
            all stages are subdirectories of project directory
        create log directory
        create main logfile

    For each stage:
        stages create/use the following files
        {stage}/{stage}.command = a list of all the executable commands for the stage
        {stage}/{stage}.complete = a list of all commands that have completed
        {stage}/{stage}.log = logfile for the stage
        {stage}/{stage}.complete = marker file, if present the stage is complete

        open/create {stage}/{stage}.log
        denovo (denovo)
            create directory {stage}
            create, save, and open command from plan
            create and open complete

        fastforward:
            if directory {stage} doesn't exist, run as for denovo
            if {stage}/finished, skip stage
            if {stage}/{stage}.command exists
                read commands
                else generate commands from plan

            if {stage}/{stage.complete exists
                remove complete commands from command and store command
                    reopen command for reading
                open complete for reading (a)


            after initial setup
                command is open for reading
                complete is open for writing (w if denovo, a if fastforward)
                stage_log is open for writing (w if denovo, a if fastforward)

        send command, complete, log to executor
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        option:         dict of options set on command line
        plan:           python object created from plan xml
        log_main        filehandle - overall logfile

        stage:          name of the current stage, from yaml
        command:        filehandle - list of commands to run in this stage
        complete:       filehandle - list of commands completed in this stage
        log_stage:      filehandle - current stage log
        -----------------------------------------------------------------------------------------"""
        self.option = getoptions()
        self.plan = {}
        self.log_main = None
        self.stage = ''
        # file pointers for commands, complete commands, and the stage log
        self.command = None
        self.complete = None
        self.log = Log()

    @staticmethod
    def open_exist(filename, mode='w'):
        """-----------------------------------------------------------------------------------------
        if filename exists, open for appending, otherwise open for writing

        :param filename: string
        :return: filehandle
        -----------------------------------------------------------------------------------------"""
        if mode == 'w':
            if os.path.exists(filename):
                mode = 'a'

        try:
            fh = open(filename, mode)
        except OSError:
            sys.stderr.write(f'Workflow:open_exist - unable to open {filename} with mode={mode}')
            exit(1)

        return fh

    def save_exist(self, filename):
        """-----------------------------------------------------------------------------------------
        if filename exists, generate filenames with a numerical suffic until an unused filename is
        found. rename the existing file to the unique name.  this allows filename to be used without
        wiping out a previous result

        :param filename: string
        :return: string             unique name for existing file
        -----------------------------------------------------------------------------------------"""
        current = self.command.name

        # add a numeric suffix to the current file and change its name
        suffix = 0
        uniquename = filename + f'.{suffix}'
        while os.path.isfile(uniquename):
            suffix += 1
            uniquename = filename + f'.{suffix}'

        os.rename(uniquename, filename)
        return uniquename

    def yaml_read(self):
        """-----------------------------------------------------------------------------------------
        read the workflow description from the yaml file

        :param filename: string     yaml workflow description
        :return: dict               workflow as a python dictionary
        -----------------------------------------------------------------------------------------"""
        fp = open(self.option['workflow'], 'r')
        self.plan = yaml.load(fp, Loader=yaml.FullLoader)
        fp.close()
        return self.plan

    def command_generate(self, stage):
        """-----------------------------------------------------------------------------------------
        generate a list of commands from the workflow for a specific stage
        1. split the command into tokens
        2. for each token, match to the other keywords in plan['stage'][stage] and replace
           %keyword with the value from the yaml. because % substitutions are processed first, they
           can include global $ definitions
        3. for each token, match to keywords in plan['definitions'] and replace $keyword with the
           value from th yaml
        4. options in <> are processed using plan['stage'][stage]['rule'] at runtime

        :param stage: string            stage in workflow
        :return: string                 command string
        -----------------------------------------------------------------------------------------"""
        current = self.plan['stage'][stage]

        # TODO check for stage finished

        command = current['command']
        token = command.split()
        # process %
        for t in range(len(token)):
            for d in current:
                if d == 'command':
                    continue
                target = f'%{d}'
                if token[t].find(target) > -1:
                    token[t] = token[t].replace(target, current[d])

        # process $
        for t in range(len(token)):
            for d in self.plan['definitions']:
                target = f'${d}'
                if token[t].find(target) > -1:
                    token[t] = token[t].replace(target, self.plan['definitions'][d])

        new = ' '.join(token)
        self.command = new
        return new

    def main_setup(self):
        """-----------------------------------------------------------------------------------------
        set up for the overall run
        1) denovo==True: create project directory, deleting current data
           denovo==False: create directory if not present, otherwise reuse current project
        2) create main log
        3) read workflow from yaml plan description

        :return:
        -----------------------------------------------------------------------------------------"""
        dir = self.option["project"]
        self.option['base'] = os.path.basename(dir)
        if self.option['denovo']:
            # in denovo mode, delete existing directory if present and create the project directory
            # and main log
            if os.path.isdir(dir):
                # project directory exists
                shutil.rmtree(dir, ignore_errors=False, onerror=None)
            os.mkdir(dir)

        else:
            # fastforward mode keeps the current directory without changing, if the project
            # directory and log file exist, proceed, otherwise create new ones

            if not os.path.isdir(dir):
                # project directory does not exist, this is a new project, create logfile below
                os.mkdir(dir)

        # start the main log
        log = f'{dir}/{self.option["base"]}.log'
        self.log.start('main', log)
        self.log.add('main', f'Begin: Project {self.option["project"]}')

        # read the workflow plan from the yaml file

        self.yaml_read()
        self.log.add('main', f'Read: Workflow {self.option["workflow"]} {len(self.plan["stage"])} stages')

        return

    def stage_setup(self, stage):
        """-----------------------------------------------------------------------------------------
        Prepare to execute commands for this stage
        1. check/create stage directory, if necessary
        2. open log file and store in object
        3. if stage is marked finished, return False, otherwise True

        if this is a denovo run, main_setup already created a new project diretory so we don't have to check


        :param dir: string      directory name
        :return: bool           False if stage is finished and should be skipped, else True
        -----------------------------------------------------------------------------------------"""
        status = True
        self.stage = stage
        self.stage_finished = False

        stagedir = f'{self.option["project"]}/{stage}'
        if not os.path.isdir(stagedir):
            # create directory if absent
            os.mkdir(dir)

        # open stage log
        stage_log = f'{w.option["project"]}/{stage}/{stage}.log'
        self.log.start('stage', stage_log)
        self.log.add('main', f'Starting stage: {self.stage} ')
        self.log.add('stage', f'Starting stage: {self.stage} ')

        # when a stage is finished a marker file {stage}/{stage}.finished is created, if this file is present
        # the stage will be skipped, return false without opening files: log.
        finished = f'{w.option["project"]}/{stage}/{stage}.finished'
        if os.path.isfile(finished):
            self.log.add(f'stage:{self.stage}: finish marker detected\n')
            self.stage_finished = True
            status = False

        return status

    def stage_finish(self):
        """-----------------------------------------------------------------------------------------
        clean up at the end of a stage
        1) mark as finished by creating stage, {base}/{stage}/{stage}.complete
        2) report stage finished in log_main
        3) close stage_log, command, and complete files


        :return:
        -----------------------------------------------------------------------------------------"""
        # mark stage finished
        finished = f'{w.option["base"]}/{stage}/{stage}.finished'
        marker = open(finished, 'r')
        marker.close()
        self.log.add('main', f'{self.stage}:finished\n')
        self.log.add('stage', f'{self.stage}:finished\n')

        for f in (self.log_stage, self.command, self.complete):
            f.close()

        return

    def stage_fast_forward(self):
        """-----------------------------------------------------------------------------------------
        When stage_fast_forward() returns False it means that commands must be generated and stored
        in the command file (done by command_generate)

        1. check if this is a --denovo run, if yes return False
        3. check if command and complete files are available, if not return False
        4. examine the list of commands {stage}.commands and completed commands {stage}.complete and
           create a new list of commands that still need to be run, return True

        :return: bool   True if commands need to be generated (denovo())
        -----------------------------------------------------------------------------------------"""
        if self.option['denovo']:
            # de novo run can't fast forward
            return False

        # try to open and read completed commands, if file is absent
        completefile = f'{w.option["base"]}/{stage}/{stage}.complete'
        done = []
        if os.path.exists(completefile):
            # read the completed job list; when run on multiple processors the complete list may not
            # be in the same order as the job list
            self.complete = open(completefile, 'r')
            for line in self.complete:
                done.append(line.rstrip())

            self.complete.close()

        # generate commands and store in commands file. for a de novo run, all files were deleted
        # so the command file will not be present

        commandfile = f'{w.option["base"]}/{stage}/{stage}.command'
        todo = []
        if os.path.exists(commandfile):
            if done:
                # since done has entries we need to fastforward
                command = open(commandfile, 'r')
                for line in self.command:
                    line = line.rstrip()
                    if line in done:
                        continue
                    todo.append(line)

                self.command.close()
            else:
                # nothing to do, the commands are generated, and none have been done
                pass
        else:
            # generate commands for this stage and save
            todo = generate_commands

        # commands have been generated and/or fastforwarded
        # if todo list is empty, mark stage finished and return False
        if not todo:
            self.stage_finish()
            return False


        # now open the new file and write the command list
        self.save_exist(commandfile)
        self.command = open(commandfile, 'w')
        for c in todo:
            self.command.write(f'{c}\n')

        # close the files so they can be reopened in read mode to run the commands
        self.command.close()
        self.complete.close()

        return True


####################################################################################################
# end of class Workflow
####################################################################################################

####################################################################################################
# command executor, runs a set of command lines in a directory,
# executor is the main part of the original rna_manager
####################################################################################################
class executor():
    """=============================================================================================
    command executor, runs a set of command lines in a directory,
    ============================================================================================="""

    def __init__(self, jobs=20):
        """-----------------------------------------------------------------------------------------
        commandlist     fh, list of command lines to run
        completelist    fh, list of completed commands
        log_stage       fh, log for details
        log_main        fh, log for general information

        jobs            number of jhob to run concurrently
        running         number of jobs currently running
        TODO maybe the main log can be done before/after running manager
        -----------------------------------------------------------------------------------------"""
        self.jobs = jobs
        self.jobid = 0
        # may not need these with new setup
        self.running = 0
        # self.total = 0
        # self.started = 0
        # self.finished = 0
        # self.succeeded = 0
        # self.failed = 0
        self.delay = 5
        self.joblist = []

    def manager(self):
        """-----------------------------------------------------------------------------------------
        main loop for executing commands in parallel.
        1. read in commands from list of commands
        2. report to log
        3. loop over manager_startjobs and manager_polljobs to run all jobs in commandlist
        TODO open file for raw error messages from subprocess
        TODO how do i know what directory.log to, probably use log directory
        -----------------------------------------------------------------------------------------"""
        # read in commands
        self.commands = []
        for line in self.commandlist:
            self.commands.append(line.rstrip())
        self.commandlist.close()

        total = len(self.commands)
        self.add('manager', 'files', f'{self.stage}',
                 f'{total} files in {self.source}')

        while self.manager_startjobs() or self.joblist:
            self.manager_polljobs()

        return self.finished

    def manager_startjobs(self):
        """-----------------------------------------------------------------------------------------
        process all commands in self.commandlist. All commands should already be complete so no
        expansion of wildcards or globs is needed

        run up to self.jobs at a time

        usage:
            while self.manager_startjobs(filelist, stage):
                self.manager_polljobs(stage)

        :param filelist: list - list of input files to process
        :param stage: dict - information describing current pipeline stage
        :return: boolean - True indicates there are more files to process
        -----------------------------------------------------------------------------------------"""
        while self.commands and self.running < self.jobs:
            # check number of jobs running and start more if necessary
            # completed commands have already been removed by fast forward, so you don't
            # have to check
            thiscommand = self.commands.pop()
            self.jobid += 1

            print(f'starting job:{self.jobid}; {thiscommand}')
            self.logstage.write('{self.stage}\t{thiscommand}')

            job = sub.Popen(thiscommand, shell=True, stdout=self.errorlog, stderr=self.errorlog)
            # TODO what is file?
            self.joblist.append([self.jobid, job, file])
            self.running += 1

        if self.commands:
            # files remain to be processed
            return True
        else:
            return False

    def manager_polljobs(self):
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
                    self.add('manager', 'complete', stage['stage'], f'jobid:{id};file:{file}')
                    self.succeeded += 1

                else:
                    # error
                    print(f'job{id} failed {file}')
                    self.add('manager', 'fail', self.current, f'status:{result};file:{file}')
                    self.failed += 1

                # include the result in the remove list, it can't be removed here because it
                # changes the list (self.joblist) that is iterating
                to_remove.append(j)

        # remove all finished jobs. Couldn't do it above because it shortens the joblist and
        # some jobs don't get polled
        for j in to_remove:
            self.joblist.remove(j)


####################################################################################################
# end of class executor
####################################################################################################

class Log(dict):
    """=============================================================================================
    Make standardized entries in log files. Only a class so it can be shared by executor and
    workflow
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        multiple logs can be created, each is associated with a filehandle by an entry in the log
        dict

        Log     dict, keys string with symbolic name for log, values file handle
        -----------------------------------------------------------------------------------------"""

    def start(self, logname, filename):
        """-----------------------------------------------------------------------------------------
        add a log to the list of logs with the key logname.  If filename exists, open the existing
        file in append mode, otherwise open a new file in write mode. If logname exists in the dict
        it is replaced.

        :param logname: string      name for log, e.g. 'main'
        :param filename: string     path to logfile
        :return: fp
        -----------------------------------------------------------------------------------------"""
        self[logname] = Workflow.open_exist(filename, 'w')
        return self[logname]

    def logtime(self):
        """-----------------------------------------------------------------------------------------
        Create a time string for use in logs year month day hour min sec concatenated

        :return: str - YYYYMoDyHrMnSc
        -----------------------------------------------------------------------------------------"""
        return time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))

    def add(self, log, message):
        """-----------------------------------------------------------------------------------------
        add a timestamped message to the named log

        :param log: string - name of log (e.g., manager, error)
        :param message: string - message
        :return: string - message written to log
        -----------------------------------------------------------------------------------------"""
        try:
            fp = self[log]
        except IndexError:
            sys.stderr.write(f'Log:add - unknown_logfile {log}\n')
            return

        # log_message = f'{tag}\t{stage}\t{message}\t{self.logtime()}\n'
        log_message = f'{self.logtime()}\t{message}\n'
        fp.write(log_message)
        fp.flush()

        return log_message


####################################################################################################
# end of class Log
####################################################################################################
####################################################################################################
# Main Program
####################################################################################################
if __name__ == '__main__':

    w = Workflow()

    now = time.localtime()
    sys.stdout.write(f'manager.py {time.asctime(now)}\n\n')
    sys.stdout.write(f'project: {w.option["project"]}\n')
    w.main_setup()
    sys.stdout.write(f'Stages read from {w.option["workflow"]}:\n')
    for stage in w.plan['stage']:
        sys.stdout.write(f'\t{stage}\n')

    for stage in w.plan['stage']:
        # create a list of commands to execute, and a file to store the list of completed commands
        w.stage_setup(stage)
        if not w.stage_fast_forward():
            # fast_forward returns false for denovo mode, or in fastforward mode if the stage
            # directory, and stage log, command, and complete files do not exist
            w.stage_denovo()

        # command execute, check for stage finished
        # TODO execute here

        w.stage_finish()

    exit(0)
