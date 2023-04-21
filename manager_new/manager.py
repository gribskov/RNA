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
        plan:           python object created from plan xml
        stage:          the current stage
        command:        the command constructed from the stage yaml
        commandlist:    filehandle - list of commands to run in this stage
        completelist:   filehandle - list of commands completeed in this stage
        option:         dict of options set on command line
        -----------------------------------------------------------------------------------------"""
        self.plan = {}
        self.stage = '
        self.command = ''
        self.commandlist = None
        self.completelist = None
        self.option = getoptions()

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

        # TODO open log, command, and complete
        logfile = f'{w.option["base"]}/{stage}/{stage}.log'
        commandfile = f'{w.option["base"]}/{stage}/{stage}.commands'
        completefile = f'{w.option["base"]}/{stage}/{stage}.complete'

        return True

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
        complete = open(completefile, 'r')
        for line in complete:
            done.append(line.rstrip())
        complete.close()

        # compare to the command list and create a new list of commands
        todo = []
        command = open(commandfile, 'r')
        for line in command:
            line = line.rstrip()
            if line in done:
                continue
            todo.append(line)
        command.close()

        # if todo_ is empty, this stage is finished, create the finished marker file and return False
        if os.path.isfile(finished):
            marker = open(finished, 'r')
            marker.close()
            return False

        # save the current command file to a new name, checking to be sure that old versions are
        # not deleted
        commandfile_old = commandfile
        suffix = 0
        commandfile_old = commandfile + f'.{suffix}'
        while os.path.isfile(commandfile_old):
            suffix += 1
            commandfile_old = commandfile + f'.{suffix}'
        os.rename(commandfile, commandfile_old)

        # now write the command list
        command = open(commandfile, 'w')
        for c in todo:
            command.write(f'{c}\n')
        command.close()

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
        w.stage_setup_dirs()
        if w.stage_fast_forward():
            # execute commands
            pass

    commands = w.command_generate('xios')
    exit(0)
