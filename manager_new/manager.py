import os.path
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

    r, restart      remove previous results and perform all stages
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

    commandline.add_argument('-r', '--restart',
                             help='Erase direcories/files and restart (%(default)s)',
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
    python manager.py <project>

    project is the base name for the run.
    For each stage is created. Each stage directory contains:
        {stage}.command = a list of all the executable commands for the stage
        {stage}.complete = a list of all commands that have completed
        {stage}.log = logfile for the stage
        {stage}.complete = marker file, if present the stage is complete

    Manager runs in two modes, controlled by --restart:
        restart: delete any existing directories/files and start a new project from scratch
        fastforward (default) use existing directories and files. Execute only stages and commands
            that haven't completed
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        option:         dict of options set on command line
        plan:           python object created from plan xml
        log_main        filehandle - overall logfile

        stage:          name of the current stage, from yaml
        command:        filehandle - list of commands to run in this stage
        complete:       filehandle - list of commands completeed in this stage
        log_stage:      filehandle - current stage log
        -----------------------------------------------------------------------------------------"""
        self.option = getoptions()
        self.plan = {}
        self.log_main = None
        self.stage = ''
        # file pointers for commands, complete commands, and the stage log
        self.command = None
        self.complete = None
        self.log_stage = None

    def open_exist(self, filename, mode='w'):
        """-----------------------------------------------------------------------------------------
        if filename exists, open for appending, otherwise open for writing

        :param filename: string
        :return: filehandle
        -----------------------------------------------------------------------------------------"""
        if mode == 'w':
            if os.path.exists(filename):
                mode = 'a'

        return open(filename, mode)

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

    @staticmethod
    def logtime():
        """-----------------------------------------------------------------------------------------
        Create a time string for use in logs year month day hour min sec
        concatenated

        :return: str - YYYYMoDyHrMnSc
        -----------------------------------------------------------------------------------------"""
        return time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))

    def main_setup(self):
        """-----------------------------------------------------------------------------------------
        set up for the overall run
        1) restart==True: create project directory, deleting current data
           restart==False: create directory if not present, otherwise reuse current project
        2) create main log
        3) read workflow from yaml plan description

        :return:
        -----------------------------------------------------------------------------------------"""
        dir = self.option["project"]
        self.option['base'] = os.path.basename(dir)
        log = f'{dir}/{self.option["base"]}.log'
        if self.option['restart']:
            # in restart mode, delete existing directory if present and create the project directory
            # and main log
            if os.path.isdir(dir):
                # project directory exists
                shutil.rmtree(dir, ignore_errors=False, onerror=None)
            os.mkdir(dir)
            self.log_main = open(log, 'w')

        else:
            # fastforward mode keeps the current directory without changing, if the project
            # directory and log file exist, proceed, otherwise create new ones

            if not os.path.isdir(dir):
                # project directory does not exist, this is a new project, create logfile below
                os.mkdir(dir)
                self.log_main = open(log, 'w')

            if os.path.isfile(log):
                # if exists, append to mexisting log
                self.log_main = open(log, 'a')
            else:
                # create log file if it doesn't exist
                self.log_main = open(log, 'w')

        self.log_main.write(f'{self.logtime()}\tBegin: Project {self.option["project"]}\n')

        # read the workflow plan from the yaml file

        self.yaml_read()
        self.log_main.write(f'{self.logtime()}\tRead: Workflow {self.option["workflow"]} ')
        self.log_main.write(f'{len(self.plan["stage"])} stages\n')

        return

    def stage_setup(self, stage):
        """-----------------------------------------------------------------------------------------
        1. check/create stage directory
        2. open log file and store in object

        :param dir: string      directory name
        :return:
        -----------------------------------------------------------------------------------------"""
        self.stage = stage
        dir = f'{self.option["project"]}/{stage}'

        if os.path.isdir(dir):
            # directory exists
            if self.option['restart']:
                # for fastforward keep the current directory without changing, otherwise
                # remove previous directory tree and create a new one
                shutil.rmtree(dir, ignore_errors=False, onerror=None)
                os.mkdir(dir)
        else:
            # create directory if absent
            os.mkdir(dir)

        # open stage log
        logfile = f'{w.option["project"]}/{stage}/{stage}.log'
        self.log_stage = self.open_exist(logfile, 'w')

        return True

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
        self.log_main.write(f'{Workflow.logtime()}\tFinished stage:{self.stage}\n')

        for f in (self.log_stage, self.command, self.complete):
            f.close()

        return

    def stage_fast_forward(self):
        """-----------------------------------------------------------------------------------------
        When stage_fast_forward() returns False it means that commands must be generated and stored
        in the command file (done by command_generate)

        1. check if this is a --restart run, if yes return False
        2. check for a file {stage}.finished, if present this stage is complete. return False
        3. check if command and complete files are available, if not return False
        4. examine the list of commands {stage}.commands and completed commands {stage}.complete and
           create a new list of commands that still need to be run, return True

        :return: bool   False if stage is complete, True if there are commands to execute
        -----------------------------------------------------------------------------------------"""
        # TODO check for restart

        # when a stage is finished a maker file {stage}.finished is created, if this file is present
        # fast forward skips the stage
        finished = f'{w.option["project"]}/{stage}/{stage}.finished'
        if os.path.isfile(finished):
            time = Workflow.logtime()
            self.log_main.write(f'{time}\tFast_forward: Finish stage:{self.stage} ')
            self.log_main.write(f'stage finish marker detected\n')
            self.log_stage.write(
                f'{time}\tFast_forward: Finish stage:{self.stage}\n')
            self.stage_finished = True
            return False
        else:
            self.stage_finished = False

        # TODO check for command and complete files

        # read the completed job list; when run on multiple processors the complete list may not
        # be in the same order as the job list

        done = []
        for line in self.complete:
            done.append(line.rstrip())
        # complete.close()

        # compare completed commands to the command list and create a new list of commands
        todo = []
        for line in self.command:
            line = line.rstrip()
            if line in done:
                continue
            todo.append(line)

        # if todo_ is empty, this stage is finished, create the finished marker file and return False
        if os.path.isfile(finished):
            self.stage_finish()
            return False

        # save the current command file to a new name, checking to be sure that old versions are
        # not deleted. the name of the current command file comes from the open filehandle
        current = self.command.name
        commandfile = f'{w.option["base"]}/{stage}/{current}'
        commandfile = current

        # add a numeric suffix to the current file and change its name
        suffix = 0
        commandfile_old = commandfile + f'.{suffix}'
        while os.path.isfile(commandfile_old):
            suffix += 1
            commandfile_old = commandfile + f'.{suffix}'

        # have to close before renamimg
        self.command.close()
        os.rename(commandfile, commandfile_old)

        # now open the new file and write the command list
        self.command = open(commandfile, 'w')
        for c in todo:
            self.command.write(f'{c}\n')

        # close the files so they can be reopened in read mode to run the commands
        self.command.close()
        self.complete.close()

        # commands need to be processed so

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
        Runs each stage keeping self.jobs running at a time.  Delay between polling is set by
        self.delay
        TODO no more handling of stages, just lists of commands
        TODO probably this whole function can be removed
        -----------------------------------------------------------------------------------------"""
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
            self.logwrite('manager', 'files', f'{stage["stage"]}',
                          f'{total} files in {self.source}')

            # check that output directories exist
            for dir in stage['dirs']:
                Pipeline.dircheck(dir)

            while self.manager_startjobs(filelist, stage) or self.joblist:
                self.manager_polljobs(stage)

        return self.finished

    def manager_startjobs(self, filelist, stage):
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
        commandlist = self.commandlist
        for command in commandlist:
            # TODO working here

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

                self.logwrite('manager', 'start', stage['stage'],
                              f'jobid:{self.jobid}; input:{file} ')
                self.logwrite('manager', 'command', stage['stage'], ' '.join(thiscommand))
                job = sub.Popen(' '.join(thiscommand), shell=True, stdout=self.errorlog,
                                stderr=self.errorlog)
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


####################################################################################################
# end of class executor
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
        if not w.stage_fast_foward():
            # fast_forward returns false for restart mode, or in fastforward mode if the stage
            # directory, and stage log, command, and complete files do not exist
            w.command_generate()

        # command execute, check for stage finished

        w.stage_finish()

    exit(0)
