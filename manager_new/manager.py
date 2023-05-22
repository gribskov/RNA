import os.path
import subprocess as sub
import shutil
import sys
import time
import argparse
import yaml
import glob


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
                             help='Erase directories/files and restart denovo (%(default)s)',
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

    Workflow runs in two modes, controlled by --restart:
    restart: delete any existing directories/files and start a new project from scratch
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
        restart (denovo)
            create directory {stage}
            create, save, and open command from plan
            create and open complete

        fastforward:
            if directory {stage} doesn't exist, run as for restart
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
                complete is open for writing (w if restart, a if fastforward)
                stage_log is open for writing (w if restart, a if fastforward)

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
        self.yaml = None
        self.log_main = None
        self.stage = ''
        self.stage_dir = ''
        self.stage_finished = False
        # file pointers for commands, complete commands, and the stage log
        self.command = None
        self.complete = None
        self.log = Log()

    @staticmethod
    def open_exist(filename, mode='w'):
        """-----------------------------------------------------------------------------------------
        if filename exists, open for appending, otherwise open for writing

        :param filename: string
        :param mode: string          file open mode, r, w, or a
        :return: filehandle
        -----------------------------------------------------------------------------------------"""
        fh = None
        file_exists = os.path.exists(filename)
        if mode == 'w':
            if file_exists:
                mode = 'a'
            fh = open(filename, mode)

        if mode == 'r':
            try:
                fh = open(filename, mode)
            except OSError:
                # make sure the file exists by opening it for writing then reopening it for reading
                # fh = open(filename, 'w')
                # fh = open(filename, 'r')
                # return False to indicate there's nothing in the file
                fh = False

        return fh

    def save_exist(self, filename):
        """-----------------------------------------------------------------------------------------
        if filename exists, generate filenames with a numerical suffix until an unused filename is
        found. rename the existing file to the unique name.  this allows filename to be used without
        wiping out a previous result

        :param filename: string
        :return: string             unique name for existing file
        -----------------------------------------------------------------------------------------"""
        # add a numeric suffix to the current file and change its name
        suffix = 0
        uniquename = filename + f'.{suffix}'
        while os.path.isfile(uniquename):
            suffix += 1
            uniquename = filename + f'.{suffix}'

        # must close file to rename
        self.command.close()
        os.rename(filename, uniquename)
        self.command = self.open_exist(filename, 'w')

        return uniquename

    def main_setup(self):
        """-----------------------------------------------------------------------------------------
        set up for the overall run
        1) restart==True: create project directory, deleting current data
           restart==False: create directory if not present, otherwise reuse current project
        2) create main log
        3) read workflow from yaml plan description

        :return:
        -----------------------------------------------------------------------------------------"""
        projdir = self.option["project"]
        self.option['base'] = os.path.basename(projdir)
        if self.option['restart']:
            # in restart mode, delete existing directory if present and create the project directory
            # and main log
            if os.path.isdir(projdir):
                # project directory exists
                shutil.rmtree(projdir, ignore_errors=False, onerror=None)
            os.mkdir(projdir)

        else:
            # fastforward mode keeps the current directory without changing, if the project
            # directory and log file exist, proceed, otherwise create new ones

            if not os.path.isdir(projdir):
                # project directory does not exist, this is a new project, create logfile below
                os.mkdir(projdir)

        # start the main log
        log = f'{projdir}/{self.option["base"]}.log'
        self.log.start('main', log)
        self.log['main'].write('\n')
        self.log.add('main', f'Project {self.option["project"]}: started')

        # expand global symbols (definitions in yaml) and store in yaml
        self.yaml = Command(filename=self.option['workflow'])
        self.yaml.read()
        self.yaml.def_main = self.yaml.expand(self.yaml.parsed['definitions'])
        self.log.add('main', f'{self.option["project"]}: workflow {self.option["workflow"]} read '
                             f'{len(self.yaml.parsed["stage"])} stages')
        sys.stdout.write(f'Stages read from {self.option["workflow"]}:\n')

        return

    def stage_setup(self, stage):
        """-----------------------------------------------------------------------------------------
        Prepare to execute commands for this stage
        1. check/create stage directory, if necessary
        2. open log file and store in object
        3. if stage is marked finished, return False, otherwise True

        if this is a restart run, main_setup already created a new project directory, so we don't have to check

        : param stage: string       name of stage to initiate
        :return: bool               False if stage is finished and should be skipped, else True
        -----------------------------------------------------------------------------------------"""
        status = True
        self.stage = stage
        self.stage_finished = False
        stagedir = self.stage_dir = f'{self.option["project"]}/{stage}'

        if not os.path.isdir(stagedir):
            # create directory if absent
            os.mkdir(stagedir)

        # open stage log
        stage_log = f'{w.option["project"]}/{stage}/{stage}.log'
        self.log.start('stage', stage_log)
        self.log.add('main', f'Stage {self.stage}: starting')
        self.log.add('stage', f'Stage started')

        # when a stage is finished a marker file {stage}/{stage}.finished is created,
        # if this file is present, the stage will be skipped,
        # return false without opening files: log.
        finished = f'{w.option["project"]}/{stage}/{stage}.finished'
        if os.path.isfile(finished):
            self.log.add('main', f'Stage:{self.stage}: skipped (finished detected)\n')
            self.log.add('stage', f'Finished marker detected\n')
            self.stage_finished = True
            self.log['stage'].close()
            status = False

        # False indicates the stage will not be run, and no processing of command/complete will
        # occur
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
        marker = open(finished, 'w')
        marker.close()

        self.log.add('stage', f'Finished marker created')
        self.log.add('main', f'Stage: {self.stage}: finished\n')
        self.log.add('stage', f'Stage finished\n')

        # close all stage files
        for f in (self.log['stage'], self.command, self.complete):
            f.close()

        return

    def stage_fast_forward(self):
        """-----------------------------------------------------------------------------------------
        When stage_fast_forward() returns False it means that commands must be generated and stored
        in the command file (done by command_generate)

        1. check if this is a --restart run, if yes return False
        3. check if command and complete files are available, if not return False
        4. examine the list of commands {stage}.commands and completed commands {stage}.complete and
           create a new list of commands that still need to be run, return True

        :return: bool   True if commands need to be generated (restart())
        -----------------------------------------------------------------------------------------"""
        self.yaml.stage = self.stage
        self.yaml.stage_dir = self.stage_dir
        # try to open and read completed commands, if file is absent the complete list will be
        # empty
        completefile = f'{w.option["base"]}/{stage}/{stage}.complete'
        done = []
        self.complete = self.open_exist(completefile, 'r')
        if self.complete:
            # read the completed job list; when run on multiple processors the complete list may not
            # be in the same order as the job list
            for line in self.complete:
                done.append(line.rstrip())
        else:
            self.complete = open(completefile, 'w')

        self.log.add('stage', f'Completed commands: {len(done)}')
        self.complete.close()
        # complete file exists and is closed

        # try to open and read commands, skip any commands in the done list
        commandfile = f'{w.option["base"]}/{stage}/{stage}.command'
        self.command = self.open_exist(commandfile, 'r')
        todo = []
        # if the command file is not present, commands must be generated from the workflow
        new_commands = False
        if not self.command:
            new_commands = True
            commands = self.yaml.command_generate()
            self.log.add('stage', f'Commands generated: {len(commands)}')

            self.command = self.open_exist(commandfile, 'w')
            for c in commands:
                self.command.write(f'{c}\n')
            self.command.close()
            self.command = self.open_exist(commandfile, 'r')
        # else:
        # command file is open and readable
        # remove any existing commands and save

        command_n = 0
        for line in self.command:
            command_n += 1
            line = line.rstrip()
            if line in done:
                continue
            todo.append(line)

        if not new_commands:
            self.log.add('stage', f'Existing commands read: {command_n}')

        self.save_exist(commandfile)
        for c in todo:
            self.command.write(f'{c}\n')

        self.log.add('stage', f'Commands (completed removed): {len(todo)}')
        self.command.close()

        # command and complete files are closed so they can be reopened by the executor

        return True


####################################################################################################
# end of class Workflow
####################################################################################################

class Command:
    """#############################################################################################
    convert the yaml file describing the workflow to a list of executable commands
    #############################################################################################"""

    def __init__(self, filename='workflow1.yaml'):
        """-----------------------------------------------------------------------------------------
        file        path to yaml file containing the workflow
        plan        python structure version of the plan
        stage       current stage
        def_main    dict, global definitions
        def_stage   dict, definitions for current stage
        -----------------------------------------------------------------------------------------"""
        self.file = filename
        self.parsed = None
        self.stage = ''
        self.stage_dir = ''
        self.def_main = {}
        self.def_stage = {}
        self.command = ''
        self.mult = {}
        self.late = {}
        self.defs = []

    def read(self):
        """-----------------------------------------------------------------------------------------
        read the workflow description from the yaml file. File is stored in Command.file

        :return: dict               workflow as a python dictionary
        -----------------------------------------------------------------------------------------"""
        fp = open(self.file, 'r')
        self.parsed = yaml.load(fp, Loader=yaml.FullLoader)
        fp.close()
        return self.parsed

    def main(self):
        """-----------------------------------------------------------------------------------------
        Set up the main definitions for the plan
        :return:
        -----------------------------------------------------------------------------------------"""
        pass

    def expand(self, defs):
        """-----------------------------------------------------------------------------------------
        recursively substitute symbols in defs with their values.

        :param defs: dict   keys are the symbols, shown as $key in commands
        :return: dict       definition dict after substitution
        -----------------------------------------------------------------------------------------"""
        stack = []
        for d in defs:
            stack.append(d)

        while stack:
            d = stack[0]
            end = len(stack)
            stack = stack[1:end]
            for thisd in defs:
                symbol = f'${thisd}'
                defs[d] = defs[d].replace(symbol, defs[thisd])
            if defs[d].find('$') > -1:
                stack.append(d)

        return defs

    def command_generate(self):
        """-----------------------------------------------------------------------------------------
        generate a list of commands from the workflow for a specific stage
        1. split the command into tokens
        2. for each token, match to the keywords in plan['stage'][stage] and replace
           $keyword with the value from the yaml.
        3. for each token, match to keywords in plan['definitions'] and replace $keyword with the
           value from th yaml
        4. options in <> are processed using plan['stage'][stage]['rule'] at runtime

        :return: string                 command string
        -----------------------------------------------------------------------------------------"""
        current = self.parsed['stage'][stage]

        # merge stage definitions with global definitions and expand the command
        # separate the definitions into the command, simple symbols, mutltiple processing
        # symbols, and late symbols
        self.def_stage = {k: self.def_main[k] for k in self.def_main}
        for item in current:
            if item == 'command':
                self.command = current[item]
            elif current[item].find('*') > -1:
                self.mult[item] = current[item]
            elif current[item].find(':') == 0:
                self.late[item] = current[item]
            else:
                self.def_stage[item] = current[item]

        self.def_stage = self.expand(self.def_stage)
        for d in self.def_stage:
            self.command = self.command.replace(f'${d}', self.def_stage[d])

        for m in self.mult:
            for d in self.def_stage:
                self.mult[m] = self.mult[m].replace(f'${d}', self.def_stage[d])
        for litem in self.late:
            for d in self.def_stage:
                self.late[litem] = self.late[litem].replace(f'${d}', self.def_stage[d])

        return self.multiple()

    def multiple(self):
        """-----------------------------------------------------------------------------------------
        based on the entries in self.mult, that is, symbols that contain wildcards, generate
        lists of matching file names to be used in individual commands

        :return:
        -----------------------------------------------------------------------------------------"""
        filelist = {}
        for symbol in self.mult:
            filelist[symbol] = glob.glob(self.mult[symbol])
            # sys.stderr.write(f'filelist:{filelist}\tsymbol[{symbol}]:{self.mult[symbol]}\n')
            if not filelist[symbol]:
                # symbol expansion  produced an empty list
                # self.log.add('stage', f'Error: expansion of symbol ${symbol} produced and empty file list')
                sys.stderr.write(f'Error: expansion of symbol ${symbol} produced and empty file list\n')

        commandlist = []
        expand = self.command
        # sys.stderr.write(f'command:{self.command}\n')
        for m in self.multiple_gen(filelist):
            # m is a dict with a value for every multiple symbol
            # sys.stderr.write(f'm:{m}\n')
            expand = self.command
            for symbol in m:
                expand = expand.replace(f'${symbol}', m[symbol])
                for litem in self.late:
                    lcom = self.late[litem][1:]
                    for symbol in m:
                        # TODO not sure why the for symbol in m loop is doubly nested, it think its so that symbols get
                        # expanded both before and after a late function runs. may not be necessary but currently works
                        # so i'm leaving as is
                        file = os.path.basename(m[symbol])
                        # sys.stderr.write(f'base:{file}\tsymbol:{symbol}\tlcom:{lcom}\n')
                        lcom = lcom.replace(f'%{symbol}', f'{file}"')
                        lcom = f'"{lcom}'
                        # sys.stderr.write(f'lcom:{lcom}\n')
                        t = eval(lcom)
                        # sys.stderr.write(f't:{t}\tl:{l}\n\n')
                        expand = expand.replace(f'${litem}', t)

            commandlist.append(expand)

        if not commandlist:
            # in case there were no symbols in command
            commandlist.append(expand)

        return commandlist

    def multiple_gen(self, filelist):
        """-----------------------------------------------------------------------------------------
        generator to create all combinations of multiple values

        :return: dict   value for each entry in self.mult
        -----------------------------------------------------------------------------------------"""
        pos = {symbol: 0 for symbol in filelist}
        while True:
            combo = {symbol: filelist[symbol][pos[symbol]] for symbol in filelist}
            yield combo
            carry = 1
            for symbol in pos:
                pos[symbol] += carry
                carry = 0
                if pos[symbol] >= len(filelist[symbol]):
                    carry = 1
                    pos[symbol] = 0

            if carry:
                break

        return


####################################################################################################
# end of class Command
####################################################################################################

class Executor:
    """#############################################################################################
    command executor, runs a set of command lines in a directory. Executor is the main part of
    the original rna_manager
    #############################################################################################"""

    def __init__(self, commandfile='stage.command', completefile='stage.complete',
                 log=None, stage='', jobs=20, delay=5):
        """-----------------------------------------------------------------------------------------
        commandfile     file of commands to open (readable)
        completefile    file of completed commands (writable)
        command         fh, list of command lines to run
        complete        fh, list of completed commands
        log             Log object, set up in Workflow

        jobs            number of job to run concurrently
        delay           number of seconds to wait between polling
        running         number of jobs currently running
        -----------------------------------------------------------------------------------------"""
        self.commandfile = commandfile
        self.command = None
        self.commandlist = []
        self.completefile = completefile
        self.complete = None
        self.log = log
        if not log:
            self.log = Log()
        self.stage = stage
        self.jobs = jobs
        self.delay = delay
        self.jobid = 0
        self.running = 0
        self.total = 0
        self.started = 0
        self.finished = 0
        self.succeeded = 0
        self.failed = 0
        self.joblist = []

    def setup(self):
        """-----------------------------------------------------------------------------------------
        main loop for executing commands in parallel.
        1. read in commands from list of commands
        2. report to log
        3. loop over manager_startjobs and manager_polljobs to run all jobs in commandlist
        -----------------------------------------------------------------------------------------"""
        # read in commands
        self.command = []
        self.command = Workflow.open_exist(self.commandfile, 'r')
        if self.command:
            for line in self.command:
                self.commandlist.append(line.rstrip())
            self.command.close()

        total = len(self.commandlist)
        self.log.add('stage', f'Executor: {total} commands to execute for stage {self.stage}')
        self.log.start('raw_error', self.log['stage'].name + '.err')

        self.complete = Workflow.open_exist(self.completefile, 'w')

        return True

    def startjobs(self):
        """-----------------------------------------------------------------------------------------
        process all commands in self.commandlist. All commands should already be complete so no
        expansion of wildcards or globs is needed

        run up to self.jobs at a time

        usage:
            while self.manager_startjobs(filelist, stage):
                self.manager_polljobs(stage)

        self.commandlist: list   list of input files to process
        self.stage: string       name of current stage in workflow
        :return: boolean                True indicates there are more files to process
        -----------------------------------------------------------------------------------------"""
        while self.commandlist and self.running < self.jobs:
            # check number of jobs running and start more if necessary
            # completed commands have already been removed by fast forward, so you don't
            # have to check
            thiscommand = self.commandlist.pop()
            self.jobid += 1
            self.log.add('stage', f'Executor: starting {thiscommand}, job ID: {self.jobid}')

            job = sub.Popen(thiscommand, shell=True,
                            stdout=self.log['raw_error'], stderr=self.log['raw_error'])
            self.log.add('stage', f'Executor: {thiscommand}')
            # job += 1
            self.joblist.append([self.jobid, job])
            self.running += 1

        if self.commandlist or self.running:
            # files remain to be processed or are still running, make sure polling continues
            return True
        else:
            return False

    def polljobs(self):
        """-----------------------------------------------------------------------------------------
        Poll the currently running jobs, and remove completed jobs from the joblist

        :return:
        -----------------------------------------------------------------------------------------"""
        # poll all jobs in joblist
        time.sleep(self.delay)
        to_remove = []
        for j in self.joblist:
            jid, job = j
            # print(f'\tjob {jid} ...', end='')
            result = job.poll()
            if result is None:
                # None indicates job is still running
                # print(f'job {jid} still running')
                pass

            else:
                # job finished
                # print(f'job {jid} finished')
                self.running -= 1
                self.finished += 1
                if result == 0:
                    # success
                    self.log.add('stage', f'Executor: complete, stage:{self.stage}, jobid:{jid}')
                    self.succeeded += 1

                else:
                    # error
                    self.log.add('stage', f'Executor: fail, stage:{self.stage}, jobid:{jid}')
                    self.failed += 1

                # include the result in the remove list, it can't be removed here because it
                # changes the list (self.joblist) that is iterating
                to_remove.append(j)

        # remove all finished jobs. Couldn't do it above because it shortens the joblist and
        # some jobs don't get polled
        for j in to_remove:
            self.joblist.remove(j)
            self.complete.write(f'{j[1].args}\n')

        return True


####################################################################################################
# end of class Executor
####################################################################################################

class Log(dict):
    """=============================================================================================
    Make standardized entries in log files. Only a class so it can be shared by Executor and
    workflow
    ============================================================================================="""

    def __init__(self, *args, **kwargs):
        """-----------------------------------------------------------------------------------------
        multiple logs can be created, each is associated with a filehandle by an entry in the log
        dict

        Log     dict, keys string with symbolic name for log, values file handle
        -----------------------------------------------------------------------------------------"""
        super(Log, self).__init__(*args, **kwargs)

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
        TODO could make this static

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

    now = time.localtime()
    sys.stdout.write(f'manager.py {time.asctime(now)}\n\n')
    w = Workflow()
    sys.stdout.write(f'Project: {w.option["project"]}\n')
    if w.option['restart']:
        sys.stdout.write(f'Restart mode (all previous files removed)\n')
    else:
        sys.stdout.write(f'Fast forward mode, continue previous run\n')

    # main set up: create project directory, start main log file
    w.main_setup()

    for stage in w.yaml.parsed['stage']:
        sys.stdout.write(f'\t{stage}\n')

    # run each stage of the workflow plan (stage: in the yaml)
    for stage in w.yaml.parsed['stage']:
        # create a list of commands to execute, and a file to store the list of completed commands
        run_this_stage = w.stage_setup(stage)
        if run_this_stage:
            w.stage_fast_forward()
            # fast_forward returns false for restart mode, or if the command file does

            commandfile = f'{w.option["base"]}/{stage}/{stage}.command'
            completefile = f'{w.option["base"]}/{stage}/{stage}.complete'
            exec = Executor(commandfile, completefile, w.log, w.stage, jobs=w.option["jobs"],
                            delay=6)
            exec.setup()
            while exec.startjobs():
                exec.polljobs()

            w.stage_finish()

    w.log.add('main', f'Project {w.option["project"]}: finished')

    exit(0)
