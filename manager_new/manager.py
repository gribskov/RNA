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
    base: project directory, all stages will be in subdirectories

    r, restart      remove previous results and perform all stages
    c, continue     continue existing run, skipping commands that are marked as complete

    options:
    j, jobs         number of concurrent jobs to run
    l, log          directory for log files

    :return: command line namespace
    -----------------------------------------------------------------------------------------"""
    base_default = './'
    commandline = argparse.ArgumentParser(
        description='Run workflow', formatter_class=arg_formatter)

    commandline.add_argument('base', type=str,
                             help='Base directory for project (%(default)s)',
                             default=base_default)

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
                             default=f'{base_default}/log/')

    args = commandline.parse_args()

    # change object values if specified on command line
    if not args.base.endswith('/'):
        args.base += '/'

    return vars(args)  # convert namespace to dict


class Workflow:

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        option:         dict of options set on command line
        plan:           python object created from plan xml
        log_main        filehandle - overall logfile
        stage:          the current stage
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

    def open_exist(self, filename):
        """-----------------------------------------------------------------------------------------
        if filename exists, open for appending, otherwise open for writing

        :param filename: string
        :return: filehandle
        -----------------------------------------------------------------------------------------"""
        mode = 'r'
        if os.path.exist(filename):
            mode = 'a'

        return = open(filename, mode)

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

    def stage_setup(self, dir):
        """-----------------------------------------------------------------------------------------
        1. check/create stage directory
        2. open log, command, complete files and store in object

        :param dir: string      directory name
        :return:
        -----------------------------------------------------------------------------------------"""
        stage = self.stage
        dir = f'{self.option["base"]}/{stage}'

        if os.path.isdir(dir):
            # directory exists
            if not self.option['fastforward']:
                # for fastforward keep the current directory without changing, otherwise
                # remove previous directory tree and create a new one
                shutil.rmtree(dir, ignore_errors=False, onerror=None)
                os.mkdirs(dir)

        # open log, command, and complete files
        logfile = f'{w.option["base"]}/{stage}/{stage}.log'
        self.log_stage = self.open_exist(logfile)
        commandfile = f'{w.option["base"]}/{stage}/{stage}.commands'
        self.command = self.open_exist(commandfile)
        completefile = f'{w.option["base"]}/{stage}/{stage}.complete'
        self.command = self.open_exist(completefile)

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
        finished = f'{w.option["base"]}/{stage}/{stage}.complete'
        marker = open(finished, 'r')
        marker.close()
        self.log_main.write(f'{Workflow.logtime()}\tFinished stage:{self.stage}\n')

        for f in (self.log_stage, self.command, self.complete):
            f.close()

        return

    def stage_fast_forward(self):
        """-----------------------------------------------------------------------------------------
        1. check for a file {stage}.finished, if present this stage is complete. return False
        2. examine the list of commands {stage}.commands and completed commands {stage}.complete and
           create a new list of commands that still need to be run, return True

        :return: bool   False if stage is complete, True if there are commands to execute
        -----------------------------------------------------------------------------------------"""

        finished = f'{w.option["base"]}/{stage}/{stage}.complete'
        if os.path.isfile(finished):
            return False

        # read in the completed job list; when run on multiple processors the complete list may not
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
        current = os.path.basename(self.command)
        commandfile = f'{w.option["base"]}/{stage}/{current}'

        # add a numeric suffix to the current file and change its name
        suffix = 0
        commandfile_old = commandfile + f'.{suffix}'
        while os.path.isfile(commandfile_old):
            suffix += 1
            commandfile_old = commandfile + f'.{suffix}'
        os.rename(commandfile, commandfile_old)

        # now open the new file and write the command list
        self.command = open(commandfile, 'w')
        for c in todo:
            self.command.write(f'{c}\n')

        # commands need to be processed so
        return True


####################################################################################################
# end of class Workflow
####################################################################################################


####################################################################################################
# Main Program
####################################################################################################
if __name__ == '__main__':

    w = Workflow()

    now = time.localtime()
    sys.stdout.write(f'manager.py {time.asctime(now)}\n\n')
    sys.stdout.write(f'project: {w.option["base"]}\n')
    w.make_directory(w.option["base"])
    w.make_directory(w.option["log"])

    # read the plan and report on the stages
    w.yaml_read()
    sys.stdout.write(f'Stages read from {w.option["workflow"]}:\n')
    for stage in w.plan['stage']:
        sys.stdout.write(f'\t{stage}\n')

    for stage in w.plan['stage']:
        # create a list of commands to execute, and a file to store the list of completed commands
        w.stage_setup_()
        if w.stage_fast_forward():
            # execute commands
            commands = w.command_generate('xios')
            pass

        w.stage_finish()

    exit(0)
