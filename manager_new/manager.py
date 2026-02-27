import os.path
# import glob
import re
import subprocess as sub
import shutil
import sys
import time
import argparse
# from os import name

import yaml
import glob
from copy import deepcopy


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

    TODO add mechanism to set environment
    :return: command line namespace converted to dict
    -----------------------------------------------------------------------------------------"""
    project_default = './'
    commandline = argparse.ArgumentParser(
        description='Run workflow', formatter_class=arg_formatter)

    # commandline.add_argument('project', type=str,
    #                          help='project directory for project (%(default)s)',
    #                          default=project_default)
    # commandline.add_argument('workflow', type=str,
    #                          help='YAML workflow plan for project (%(default)s)',
    #                          default='workflow.yaml')

    commandline.add_argument('-r', '--restart',
                             help='Erase current directories/files and restart denovo (%(default)s)',
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
    # if not args.project.endswith('/'):
    #     args.project = args.project.rstrip('/')

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
        option:         dict of options set on command line (from argparse())
        command:        Command object, parsed commands and stage information and command templates
        log:            Log object - all logs
        -----------------------------------------------------------------------------------------"""
        self.option = None
        self.command = None
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
                # TODO error message
                fh = False

        return fh

    def save_exist(self, filename):
        """-----------------------------------------------------------------------------------------
        if filename exists, generate filenames with a numerical suffix until an unused filename is
        found. rename the existing file to the unique name.  this allows filename to be used without
        wiping out a previous result

        # TODO make static

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

    def prepare_project(self):
        """-------------------------------------------------------------------------------------------------------------
        1) read the workflow yaml file
        2) check if this is a restart run
            if  restart==True:
            create project directory, deleting current data

            elif restart==False:
            create directory if not already present, otherwise reuse current directory

        3) if necessary, create directories specified in directories: section of workflow
        4) create/open stdin, stdout, completed_commands, and log files
        5) resolve symbols in commands in commands: section of workflow
           store command templates

        :return:
        -------------------------------------------------------------------------------------------------------------"""
        # read and parse workflow yaml
        self.command = Command(filename=self.option['workflow'])
        project = self.command.parsed['project']

        # create directories
        self.dir_exist(project)
        for d in self.command.parsed['directories']:
            dir = f'{project}/{d}'
            self.dir_exist(dir)

        # create log object and add main, stderr, and stdout logs. all have the same basename as the project
        # log.start will not delete log if it already exists.
        project = self.command.parsed['project']
        log_base_name = f'{project}/{os.path.basename(project)}'
        log = self.log = Log()
        log.start('main', log_base_name + '.log')
        log.start('stdout', log_base_name + '.out')
        log.start('stderr', log_base_name + '.err')
        log.add('main', f'Project {project}: started')
        log.add('main', f'{project}: workflow {self.option["workflow"]} read '
                        f'{len(self.command.parsed["commands"])} commands')

        # check for restart mode, project directory must be created or identified before starting logs
        if self.option['restart']:
            # in restart mode, delete the existing directory tree, if present ,and create the project directory
            log.add('main', f'Restart mode. Project directory ({project}) will be replaced')
            if os.path.isdir(project):
                # project directory exists, remove the entire directory tree
                shutil.rmtree(project, ignore_errors=False, onerror=None)
            os.mkdir(project)
        elif os.path.isdir(project):
            # project directory exists and will be added to
            log.add('main', f'Project directory ({project}) exists and will be reused')
        else:
            # project directory does not exist, this is a new project
            log.add('main', f'Project directory ({project}) created')
            os.mkdir(project)

        # expand static symbols
        # setup templates

        return None

    @staticmethod
    def dir_exist(dirpath):
        """-------------------------------------------------------------------------------------------------------------
        checks if directory path exists, and if it does not, creates it

        :param dirpath: string      path to directory
        :return: None
        -------------------------------------------------------------------------------------------------------------"""
        if not os.path.isdir(dirpath):
            # create directory if absent
            os.mkdir(dirpath)

        return None

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
        # TODO change detection of stage start
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

    # def stage_finish(self):
    #     """-----------------------------------------------------------------------------------------
    #     clean up at the end of a stage
    #     1) mark as finished by creating stage, {base}/{stage}/{stage}.complete
    #     2) report stage finished in log_main
    #     3) close stage_log, command, and complete files
    #
    #     :return:
    #     -----------------------------------------------------------------------------------------"""
    #     # mark stage finished
    #     finished = f'{w.option["base"]}/{stage}/{stage}.finished'
    #     marker = open(finished, 'w')
    #     marker.close()
    #
    #     self.log.add('stage', f'Finished marker created')
    #     self.log.add('main', f'Stage: {self.stage}: finished\n')
    #     self.log.add('stage', f'Stage finished\n')
    #
    #     # close all stage files
    #     for f in (self.log['stage'], self.command, self.complete):
    #         f.close()
    #
    #     return

    def stage_fast_forward(self):
        """-----------------------------------------------------------------------------------------
        TODO defer fast_forward until after main command loop is updated
        When stage_fast_forward() returns False it means that commands must be generated and stored
        in the command file (done by command_generate)

        1. check if this is a --restart run, if yes return False
        3. check if command and complete files are available, if not return False
        4. examine the list of commands {stage}.commands and completed commands {stage}.complete and
           create a new list of commands that still need to be run, return True

        :return: bool   True if commands need to be generated (restart())
        -----------------------------------------------------------------------------------------"""
        self.command.stage = self.stage
        self.command.stage_dir = self.stage_dir
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
            commands = self.command.command_generate()
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

class Template:
    """#############################################################################################
    A template is the command for a specific stage with all static variables completed and wildcards
    converted to symbols (marked by %) that can be filled in during the run (see Template.fill())
    #############################################################################################"""

    def __init__(self, name='', priority=0, command=''):
        """-----------------------------------------------------------------------------------------
        command - expanded command with static symbols ($ symbols) expanded
        priority - integer giving a priority to this template, 0 is first
        realtime - symbols in the process of expansion and execution (formerly called "late"
                   symbols)
        status - 'not started', 'running', 'finished'
        -----------------------------------------------------------------------------------------"""
        self.name = name
        self.command = command
        self.priority = priority
        self.realtime = {}
        self.created = set()
        self.status = 'not started'
        self.dir = ''
        self.error = ''
        self.log = ''

    def fill(self):
        """-----------------------------------------------------------------------------------------
        Prepare the expanded command for execution with runtime lists of files
        for each filename with wildcard, create list of matching files

        :%in.set(str) defines a glob that matches symbol %in and str
        :%in.replace(from,to) creates a new filename from glob %in

        :return:
        -----------------------------------------------------------------------------------------"""
        # TODO make these class variables and precompile
        setre = r'(?P<expression>%(?P<symbol>[^.]+).set\((?P<value>[^)]+)\))'
        replacere = r'(?P<expression>%(?P<symbol>[^.]+).replace\((?P<old>[^,]+),\s*(?P<new>[^)]+)\))'

        command = self.command
        command_list = []
        target = []
        if command.find('%') == -1:
            # command has no % expressions to process, just add to the commandlist and return
            if command not in self.created:
                command_list.append({'command': command, 'priority': self.priority, 'stage': self.name})
                self.created.add(command)

            return command_list

        for m in re.finditer(setre, command):
            # find % set expressions, target are the files that match the glob in the set expression, e.g. *.xios
            print(f'globbing:{m.group("value")}')
            target = Template.glob_update(m.group('value'))

        for t in target:
            # replace % set expression with targets using % replace expression
            # check to see this target has not been previously  processed
            # add to processed list
            result = re.sub(setre, t, command)
            print(f'result:{result}')
            # get basename of globbed filename
            basetarget = os.path.basename(t)
            if basetarget in self.created:
                # skip globbed names that have already been processed
                continue

            # process % replace expression
            for m in re.finditer(replacere, command):
                # find % replace expressions, there may be more than one
                print(f'expression"{m.group("expression")} symbol:{m.group("symbol")} '
                      f'replacement:{m.group("old")} => {m.group("new")}')
                old = m.group('old').replace('"', '').replace("'", "")
                new = m.group('new').replace('"', '').replace("'", "")
                renamed = basetarget.replace(old, new)
                print(f'before:{basetarget}\t\tafter:{renamed}')
                result = result.replace(m.group('expression'), renamed)
                print(f'before:{basetarget}\t\tafter:{result}')
                command_list.append({'stage': self.name, 'priority': self.priority, 'command': result})

            # mark the target as processed
            self.created.add(basetarget)

        return command_list

    @staticmethod
    def glob_update(pattern):
        """-----------------------------------------------------------------------------------------
        match the glob pattern and return a list of matching filenames

        :param pattern: str     a globbing pattern corresponding to a filepath
        :return: list           matching files
        -----------------------------------------------------------------------------------------"""
        # re engine does not like \\ so convert to /
        return [m.replace('\\', '/') for m in glob.glob(pattern)]


class Command:
    """#################################################################################################################
    convert the yaml workflow file describing the workflow to a list of executable commands

    each stage can add/replace symbols in the global definitions, and can have tokens that can only
    be expanded at run time (see class Template())

    TODO add description of workflow syntax
    #################################################################################################################"""

    def __init__(self, filename='workflow1.yaml'):
        """-------------------------------------------------------------------------------------------------------------
        filename        path to yaml file containing the workflow
        parsed          yaml parsed to python dictionary
        template        list Template
                        command templates commands prepared for populating with filenames at run time
        static_symbols  symbols defined in workflow and referenced as $symbol.
                        symbols can be expanded before the run starts and do not change during the run
                        TODO not exactly true, they can differ in different commands
        -------------------------------------------------------------------------------------------------------------"""
        self.filename = filename
        self.parsed = None
        self.templates = []
        self.static_symbols = {}
        self.commands = []
        if filename:
            # expand static commands and create templates if filename exists
            self.parse_workflow()

    def make_templates(self):
        """-------------------------------------------------------------------------------------------------------------
        build the command templates. Templates have all static symbols expanded and any symbols to be
        determined at run time replaced with % symbols. template is a dict with the command names as keys, directly from
        the parsed workflow yaml

        The workflow yaml must have been read and parsed and stored in self.parsed
        :return:
        -------------------------------------------------------------------------------------------------------------"""
        priority = 0
        if not self.parsed:
            # error yaml not parsed
            pass
            # TODO complete this error to error log

        for com in self.parsed['commands']:
            # expand stage-specific symbols
            expanded_stage_defs = self.expand(com)
            template = Template(com, priority, expanded_stage_defs['$command'])
            self.templates.append(template)
            priority += 1

            return len(self.templates)

    def parse_workflow(self):
        """-------------------------------------------------------------------------------------------------------------
        read workflow from yaml file and expand global static symbols
        expand static symbols

        keep this as a separate method in case it is not convenient to parse the workflow file when creating the Command
        object.

        :return: dict   parsed yaml as python dictionary
        -------------------------------------------------------------------------------------------------------------"""
        # read the workflow yaml file
        self.read_parse()
        # expand global static symbols, do not process stage commands at this point
        self.static_symbols = self.expand(self.parsed['definitions'])

        return self.parsed

    def read_parse(self):
        """-----------------------------------------------------------------------------------------
        read the workflow description from the yaml file.

        :return: dict               workflow as a python dictionary
        -----------------------------------------------------------------------------------------"""
        fp = open(self.filename, 'r')
        try:
            self.parsed = yaml.load(fp, Loader=yaml.FullLoader)
        except yaml.YAMLError:
            pass
            # TODO complete this error to error log

        fp.close()

        return self.parsed

    def expand(self, definitions):
        """-----------------------------------------------------------------------------------------
        Expand all static symbols (indicated by $ in workflow). Return symbols with the
        $symbol as the key

        :param definitions: dict    definition lines converted from workflow yaml
        :return: dict               expanded symbols dict
        -----------------------------------------------------------------------------------------"""
        stack = []
        symbol = deepcopy(self.static_symbols)
        # Symbols may be defined in terms of other symbols so recursive expansion
        # is necessary. As translated from workflow yaml defs is a dictionary with the symbols as
        # keys (no $)
        for d in definitions:
            # add $ to the definition keys and push definitions all definitions on stack
            stack.append(['$' + d, definitions[d]])
        # TODO need to process already defined symbols first so they override global symbols
        while stack:
            # recursively expand symbols, symbols whose values do not contain $ are moved to symbol
            # TODO should check for unexpandable symols, maybe stack not changing in size for multiple rounds?
            dkey, dval = stack.pop()
            if dval.find('$') == -1:
                # no expandable symbol, save to symbol dict. this is the only place definitions are
                # removed from stack
                symbol[dkey] = dval
                continue

            for s in symbol:
                # expandable symbol detected, check all currently known symbols to see which it is
                # (or unknown)
                if dval.find(s) == -1:
                    # symbol not present, skip to next symbol
                    continue

                # symbol s is in this definition, expand it
                if any(c in dval for c in '*?'):
                    # wildcard replace with % expression, escape * and ? to prevent duplicating
                    # escape replaces c with &ord(c);
                    print(dval)
                    d = dkey.replace('$', '')
                    dval = dval.replace('*', f'&{ord("*")};')
                    dval = dval.replace('?', f'&{ord("?")};')
                    dval = f'%{d}.set({dval})'
                else:
                    dval = dval.replace(s, symbol[s])

            # if you reach here the symbol is either unknown (needs more expansion) or has been
            # completely expanded (next time it is popped it will go into the symbol dict)
            stack = [[dkey, dval]] + stack

        # convert escaped symbols back to characters
        r1 = f'&{ord("*")};'
        r2 = f'&{ord("?")};'
        for s in symbol:
            symbol[s] = symbol[s].replace(r1, '*').replace(r2, '?')

        return symbol

    def generate(self):
        """-----------------------------------------------------------------------------------------
        generate executable commands from stage templates

        updating to generate commands for all stages

        generate a list of commands from the workflow for a specific stage
        1. split the command into tokens
        2. for each token, match to the keywords in plan['stage'][stage] and replace
           $keyword with the value from the yaml.
        3. for each token, match to keywords in plan['definitions'] and replace $keyword with the
           value from th yaml
        4. options in <> are processed using plan['stage'][stage]['rule'] at runtime

        break command lines from workflow yaml into
        commands
        mult - symbols generated with wild cards
        late - dependent symbols created from other symbols. they are late because they generally
                can't be resolved until the previous stage finishes

        :return: string                 command string
        -----------------------------------------------------------------------------------------"""
        for stage in self.templates:
            # stage is a Template object
            print(stage.name)
            # stage_symbol = self.expand(self.parsed['stage'][stagename])
            # all symbols have been expanded so the only thing we need is the final command for the stage
            # this_stage = Template(name=stagename, command=stage_symbol['$command'])
            self.commands += stage.fill()
            # print(stage)

        return len(self.commands)

    def multiple(self):
        """-----------------------------------------------------------------------------------------
        based on the entries in self.mult, that is, symbols that contain wildcards, generate
        lists of matching file names to be used in individual commands

        TODO review, may no longer be used

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
        TODO review, may no longer be used

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
    command executor, runs a set of command lines from Command.commands. each item in the command
    list is created in command.generate() using Template.fill()
    {'stage': stage.name, 'priority': priority, 'command': full command}
    #############################################################################################"""

    def __init__(self, commandlist=None, log=None, jobs=20, delay=5):
        """-----------------------------------------------------------------------------------------
        commandlist     Command.commands - list of commands to execute
        # completefile    file of completed commands (writable)
        command         fh, list of command lines to run
        complete        fh, list of completed commands
        log             Log object, set up in Workflow

        jobs            number of job to run concurrently
        delay           number of seconds to wait between polling
        running         number of jobs currently running
        -----------------------------------------------------------------------------------------"""
        # self.commandfile = commandfile
        # self.command = None
        self.commandlist = commandlist
        # self.completefile = completefile
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
        # the complete commands are now stored in a log called projectname.complete
        # self.complete = Workflow.open_exist(self.completefile, 'w')

        return True

    def prioritize(self):
        """-----------------------------------------------------------------------------------------
        sort command list by stage priority
        :return:
        -----------------------------------------------------------------------------------------"""
        commandlist = self.commandlist.commands
        self.commandlist.commands = sorted(commandlist, key=lambda c: c['priority'])

        return

    def startjobs(self):
        """-----------------------------------------------------------------------------------------
        start self.jobs commands in self.commandlist. All commands should already be complete so no
        expansion of wildcards or globs is needed

        self.commandlist: list   list of input commands to process
        :return: boolean         True indicates there are more files to process
        -----------------------------------------------------------------------------------------"""
        commandlist = self.commandlist.commands
        while commandlist and self.running < self.jobs:
            # check number of jobs running and start more if necessary
            # completed commands have already been removed by fast forward, so you don't
            # have to check
            entry = commandlist[0]
            commandlist = commandlist[1:]
            thiscommand = entry['command']
            self.jobid += 1
            self.log.add(entry['stage'] + 'log', f'Executor: starting {thiscommand}, job ID: {self.jobid}')
            job = sub.Popen(thiscommand, shell=True,
                            env={'DATAPATH': '/scratch/scholar/mgribsko/RNAstructure/data_tables'},
                            stdout=self.log['stdout'], stderr=self.log['stderr'])
            self.log.add(entry['stage'] + 'log', f'Executor: {thiscommand}')
            # job += 1
            self.joblist.append([self.jobid, job, entry['stage']])
            self.running += 1

        if self.joblist or self.running:
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
            jid, job, stage = j
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
                    self.log.add(stage + 'log', f'Executor: complete, jobid:{jid} stage:{stage}, ')
                    self.succeeded += 1

                else:
                    # error
                    self.log.add(stage + 'err', f'Executor: fail, jobid:{jid} stage:{stage}, ')
                    self.failed += 1

                # include the result in the remove list, it can't be removed here because it
                # changes the list (self.joblist) that is iterating
                to_remove.append(j)

        # remove all finished jobs. Couldn't do it above because it shortens the joblist and
        # some jobs don't get polled
        for j in to_remove:
            self.joblist.remove(j)
            message = f'status={j[1].returncode}\t{j[1].args}'
            if j[1].returncode == 0:
                self.log.add('complete', message)
            else:
                # returned error status != 0
                self.log.add('main', 'failed' + '\t' + message)

            # self.complete.write(f'{j[1].args}\n')

        return True


####################################################################################################
# end of class Executor
####################################################################################################

class Log(dict):
    """=============================================================================================
    Make standardized entries in log files. A class so it can be shared by Workflow, Command, and
    Executor
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
        return time.strftime('%Y.%m.%d.%H%M%S', time.localtime(time.time()))

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
            log_message = f'Log:add - unknown_logfile {log}\n'
            sys.stderr.write(log_message)
            return log_message

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
    w.option = getoptions()
    # main set up: read workflow, create project directory, start main log file
    w.prepare_project()

    # sys.stdout.write(f'Project: {w.project}\n')
    # if w.option['restart']:
    #     sys.stdout.write(f'Restart mode (all previous files removed)\n')
    #     w.log.add('main', f'Project {w.project}: restart, directory reset')
    # else:
    #     sys.stdout.write(f'Fast forward mode, continue previous run\n')
    #     w.log.add('main', f'Project {w.project}: fastforward, skipping completed commands')

    #
    # main loop over all complete commands
    #
    w.command.generate()
    exec = Executor(w.command, w.log, jobs=w.option["jobs"], delay=4)
    while exec.commandlist.commands:
        # continue as long as there are commands in the command list
        exec.prioritize()
        exec.startjobs()
        exec.polljobs()

        # generate more commands if possible
        w.command.generate()

    # run each stage of the workflow plan (stage: in the yaml)
    # for stage in w.yaml.parsed['stage']:
    #     # create a list of commands to execute, and a file to store the list of completed commands
    #     run_this_stage = w.stage_setup(stage)
    #     if run_this_stage:
    #         w.stage_fast_forward()
    #         # fast_forward returns false for restart mode, or if the command file does
    #
    #         commandfile = f'{w.option["base"]}/{stage}/{stage}.command'
    #         completefile = f'{w.option["base"]}/{stage}/{stage}.complete'
    #         exec = Executor(commandfile, completefile, w.log, w.stage, jobs=w.option["jobs"],
    #                         delay=4)
    #         exec.setup()
    #         # TODO job execution turned off
    #         # while exec.startjobs():
    #         #     exec.polljobs()
    #
    #         w.stage_finish()

    w.log.add('main', f'Project {w.option["project"]}: finished')

    exit(0)
