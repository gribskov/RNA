import argparse
import glob
import os.path
import re
import shutil
import subprocess
import sys
import time
from copy import deepcopy
import yaml


def arg_formatter(prog):
    """---------------------------------------------------------------------------------------------
    Set up formatting for help. Used in arg_get.

    :param prog:
    :return: argparse formatter class
    ---------------------------------------------------------------------------------------------"""
    return argparse.HelpFormatter(prog, max_help_position=60, width=120)


def getoptions():
    """---------------------------------------------------------------------------------------------
    Command line arguments:
        required positional argument: workflow file

    options:
    j, jobs         number of concurrent jobs to run
    l, log          directory for log files
    q, quiet        minimal output to terminal
    r, restart      remove previous results and perform all stages

    TODO add mechanism to set environment variables
    :return: command line namespace converted to dict
    ---------------------------------------------------------------------------------------------"""
    project_default = './'
    commandline = argparse.ArgumentParser(
        description='Run workflow', formatter_class=arg_formatter)

    # commandline.add_argument('project', type=str,
    #                          help='project directory for project (%(default)s)',
    #                          default=project_default)
    commandline.add_argument('workflow', type=str,
                             help='YAML workflow plan for project (%(default)s)',
                             default='workflow.yaml')

    commandline.add_argument('-r', '--restart',
                             help='Erase current directories/files and restart denovo (%(default)s)',
                             default=False, action='store_true')

    # commandline.add_argument('-w', '--workflow', type=str,
    #                          help='YAML workflow plan for project (%(default)s)',
    #                          default=f'workflow.yaml')

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
        self.project = ''
        # print(f'python version: {sys.version}')

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
                # this error likely occurs when opening log files so write error to STDERR
                sys.stderr.write(f'Workflow:open_exist - unable to open ({filename})')
                fh = False

        return fh

    def save_exist(self, filename):
        """-----------------------------------------------------------------------------------------
        if filename exists, generate filenames with a numerical suffix until an unused filename is
        found. rename the existing file to the unique name.  this allows filename to be used without
        wiping out a previous result

        :param filename: string     path to file
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
        """-----------------------------------------------------------------------------------------
        1) read the workflow yaml file and expand symbols
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
        -----------------------------------------------------------------------------------------"""
        # read and parse workflow yaml, expand static symbols
        self.command = Command(filename=self.option['workflow'])
        project = self.project = self.command.parsed['project']
        use_fast_forward = True

        # check for restart mode, project directory must be created/identified before starting logs
        restart_comment = f'Project directory ({project}) created'
        if self.option['restart']:
            # in restart mode, delete the existing directory tree, if present
            # log.add('main',
            use_fast_forward = False
            if os.path.isdir(project):
                # project directory exists, remove the entire directory tree
                shutil.rmtree(project, ignore_errors=False, onerror=None)
                restart_comment = f'Project {project}\tRestart mode. Project directory ({project}) will be replaced'

        elif os.path.isdir(project):
            # directory already exists and will be used
            restart_comment = f'Project {project}\tProject directory ({project}) exists and will be reused'

        # create directories, dir_exist will not create a directory if it already  exists
        self.dir_exist(project)
        for d in self.command.parsed['directories']:
            thisdir = f'{project}/{d}'
            self.dir_exist(thisdir)

        # create log object and add main, stderr, and stdout logs. all have the same basename as
        # the project
        # log.start will not delete log if it already exists.
        log_base_name = f'{project}/{os.path.basename(project)}'
        log = self.log = Log()
        self.command.log = log
        log.start('main', log_base_name + '.log')
        log.start('stdout', log_base_name + '.out')
        log.start('stderr', log_base_name + '.err')
        log.add('main', f'Project {project}: started')
        log.add('main', restart_comment)
        log.add('main', f'Project {project}\tworkflow\tread {len(self.command.parsed["commands"])} commands')

        # setup templates
        template_n = self.command.make_templates()

        # fast forward
        if use_fast_forward:
            self.fast_forward()

        return template_n

    @staticmethod
    def dir_exist(dirpath):
        """-----------------------------------------------------------------------------------------
        checks if directory path exists, and if it does not, creates it

        :param dirpath: string      path to directory
        :return: None
        -----------------------------------------------------------------------------------------"""
        if not os.path.isdir(dirpath):
            # create directory if absent
            os.mkdir(dirpath)

        return None

    def fast_forward(self):
        """-----------------------------------------------------------------------------------------
        If using an existing project, fast_forward skips previously completed commands and begins
        at the following commands.
        1) look up started and finished commands in project.log
        2) generate the template.created lists for the completed commands (all started commands)
        3) add commands that were started but not finished to the current list

        :return:
        -----------------------------------------------------------------------------------------"""
        mainlog = self.log['main'].name
        main = open(mainlog, 'r')

        started = []
        finished = []
        for line in main:
            if 'Executor\tstarted' in line:
                started.append(line.rstrip())
            if 'Executor\tfinished' in line:
                # need the command to check against jobs started
                logtime, msgsource, msg, template, command, jobid = line.rstrip().split('\t')
                finished.append(command)

        # make a dictionary of the created list for each template
        tidx = {}
        for t in self.command.templates:
            tidx[t.name] = t

        command_list = []
        for command in started:
            logtime, msgsource, msg, template, command, jobid = command.split('\t')
            # print(time, msgsource, msg, template, command, jobid)
            tidx[template].created.add(command)
            if command not in finished:
                # the command list format has to follow Template.fill()
                command_list.append({'command': command,
                                     'priority': tidx[template].priority,
                                     'commandname': template})

        self.command.commands = command_list
        return command_list


####################################################################################################
# end of class Workflow
####################################################################################################

class Template:
    """#############################################################################################
    A template is the command for a specific stage with all static variables completed and wildcards
    converted to symbols (marked by %) that can be filled in during the run (see Template.fill())
    #############################################################################################"""
    # regular expressions for parsing workflow
    setre = re.compile(r'(?P<expression>%(?P<symbol>[^.]+).set\((?P<value>[^)]+)\))')
    replacere = re.compile(r'(?P<expression>%(?P<symbol>[^.]+).replace\((?P<old>[^,]+),\s*(?P<new>[^)]+)\))')

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
        # self.dir = ''
        # self.error = ''
        # self.log = ''

    def fill(self):
        """-----------------------------------------------------------------------------------------
        Prepare the expanded command for execution with runtime lists of files
        for each filename with wildcard, create list of matching files

        %in.set(pathglob) defines a glob that replaces symbol %in and glob(pathglob), e.g.,
            %in.set(/scratch/bell/mgribsko/rna/data/curated/fasta/rnasep*.fa)
        %in.replace(from,to) creates a new filename from glob %in, e.g.,
            %in.replace('.fa', '.pfs')

        :return: list       completed commands (ready to execute)
        -----------------------------------------------------------------------------------------"""
        command = self.command
        command_list = []
        target = []
        if command.find('%') == -1:
            # command has no % expressions to process, just add to the command_list and return
            if command not in self.created:
                command_list.append({'command': command,
                                     'priority': self.priority,
                                     'stage': self.name})
                self.created.add(command)

            return command_list

        for m in Template.setre.finditer(command):
            # find % set expressions, targets are the files that match the glob in the
            # set expression, e.g. *.xios
            target = Template.glob_update(m.group('value'))

        for t in target:
            # replace % set expression with targets
            result = Template.setre.sub(t, command)
            # get basename of globbed filename
            basetarget = os.path.basename(t)

            # process % replace expression
            for m in Template.replacere.finditer(command):
                # find % replace expressions, there may be more than one
                # print(f'expression"{m.group("expression")} symbol:{m.group("symbol")} '
                #       f'replacement:{m.group("old")} => {m.group("new")}')
                # remove quotes in re matches
                old = m.group('old').replace('"', '').replace("'", "")
                new = m.group('new').replace('"', '').replace("'", "")
                renamed = basetarget.replace(old, new)
                result = result.replace(m.group('expression'), renamed)
                if result not in self.created:
                    # prevent completing commands more than once
                    command_list.append({'commandname': self.name,
                                         'priority': self.priority,
                                         'command': result})
                    self.created.add(result)

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
    """#############################################################################################
    convert the yaml workflow file describing the workflow to a list of executable commands

    each stage can add/replace symbols in the global definitions, and can have tokens that can only
    be expanded at run time (see class Template())

    1) project: Project name (required) used as the name for the output directory
       use -r to completely remove all existing content
    2) definitions: Symbols used in directories: and commands sections
       these can be paths or aliases for executables
    3) directories: Directories to create under the project directory
       usually used for output
    4) commands: Named template commands to be run, in order of listing
          command: (required)
          wildcards are allowed as inputs
          %replace(old,new) generates output filenames from inputs

   example:
    # workflow for calculating fingerprints from sequence using partition function
   ---
  project: testmanager
  definitions:
    #   code locations
    python: python
    RNAstructure: ../../RNAstructure/exe
    XIOSexe: ..
    fasta: /scratch/scholar/mgribsko/standards/fasta_fixed
  directories:
    # directories to create under project; also creates symbols
    partition: $project/partition
    stochastic: $project/stochastic
    xios: $project/xios
  commands:
    partition:
      command: $RNAstructure/partition $in $out $option
      option:
        -t 280
        -q
      in: $fasta/rnasep_ar*.fa
      out: $partition/%in.replace('.fa', '.pfs')
    stochastic:
      command: $RNAstructure/stochastic $in $out $option
      option: -s 3
      in: $partition/*.pfs
      out: $stochastic/%in.replace('.pfs', '.ct')
    xios:
      command: $python $XIOSexe/stochastic_to_xios.py $option $in $out
      option: -m 3 -c 50
      in: $stochastic/*.ct
      out: $xios/%in.replace('.ct', '.xios')
    #############################################################################################"""

    def __init__(self, filename='workflow1.yaml'):
        """-----------------------------------------------------------------------------------------
        filename        path to yaml file containing the workflow
        parsed          yaml parsed to python dictionary
        template        list Template
                        command templates prepared for populating with filenames at run time
        static_symbols  symbols defined in workflow and referenced as $symbol.
                        symbols are expanded before the run starts and do not change during the run,
                           unless overidden in one of the named commands
        -----------------------------------------------------------------------------------------"""
        self.filename = filename
        self.parsed = None
        self.templates = []
        self.static_symbols = {}
        self.commands = []
        self.log = None
        if filename:
            # expand static commands and create templates if filename exists
            self.parse_workflow()

    def make_templates(self):
        """-----------------------------------------------------------------------------------------
        build the command templates. Templates have all static symbols expanded and any symbols to
        be determined at run time replaced with % symbols. template is a dict with the command names
        as keys, directly from the parsed workflow yaml

        The workflow yaml must have been read and parsed and stored in self.parsed

        :return: int    number of templates read
        -----------------------------------------------------------------------------------------"""
        priority = 0
        if not self.parsed:
            # error yaml not parsed
            self.log.add('main', f'Project\tError\tparsed yaml not found\t\t')

        commands = self.parsed['commands']
        for com in commands:
            # symbols can occur in the command section, expand these too
            local_defs = self.expand(commands[com])
            template = Template(com, priority, local_defs['$command'])
            self.templates.append(template)
            # TODO allow commands to override default priority
            priority += 1

        return len(self.templates)

    def parse_workflow(self):
        """-----------------------------------------------------------------------------------------
        read workflow from yaml file and expand global static symbols
        expand static symbols

        keep this as a separate method in case it is not convenient to parse the workflow file when
        creating the Command object.

        :return: dict   parsed yaml as python dictionary
        -----------------------------------------------------------------------------------------"""
        # read the workflow yaml file
        self.read_parse()
        # expand global static symbols, do not process commands at this point
        # process directories, then definitions
        self.static_symbols = {'$project': self.parsed['project']}
        self.static_symbols.update(self.expand(self.parsed['directories']))
        self.static_symbols.update(self.expand(self.parsed['definitions']))

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
            self.log.add('main', f'Project\tError\tyaml not loaded\t\t')

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

        while stack:
            # recursively expand symbols, symbols whose values do not contain $ are moved to symbol
            # TODO should check for unexpandable symbols, maybe stack not changing in size for
            #  multiple rounds?
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
                    # print(dval)
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
        generate executable commands from command templates
        for each command template in Command.templates, complete the template by adding in files
        that are available only at runtime. These are indicated by %variable_name

        :return: string                 command string
        -----------------------------------------------------------------------------------------"""
        for template in self.templates:
            # template is a Template object
            newcommands = template.fill()
            self.commands += newcommands
            if newcommands:
                # only log if new commands are generated
                self.log.add('main', f'Command\tGenerate commands\t{len(newcommands)} '
                                     f'commands added\t{template.name}\t')

        return len(self.commands)


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
        self.commandlist = commandlist
        self.joblist = []
        self.complete = None
        self.log = log
        if not log:
            self.log = Log()
        self.jobs = jobs
        self.delay = delay
        self.jobid = 0
        self.finished = 0
        self.succeeded = 0
        self.running = 0
        # self.total = 0
        # self.started = 0
        self.failed = 0

    def prioritize(self):
        """-----------------------------------------------------------------------------------------
        sort command list by stage priority
        :return:
        -----------------------------------------------------------------------------------------"""
        commandlist = self.commandlist.commands
        self.commandlist.commands = sorted(commandlist, key=lambda c: c['priority'], reverse=True)

        return

    def conda_run(self, cmd):
        """-----------------------------------------------------------------------------------------
        Run a command as a subprocess in a conda environment
        For a python script, the command could be 'python your_script.py'

        :param cmd: string      command to run
        :return: int            job id
        -----------------------------------------------------------------------------------------"""
        # Use the CONDA_PREFIX environment variable to target the current environment
        # TODO this works but it needs to not have the data tables hardwired
        env = os.environ.copy()
        env["DATAPATH"] = ('/scratch/bell/mgribsko/rna/RNAstructure/data_tables')
        project_dir = os.getenv("PROJECT_DIR")
        conda_prefix = os.environ.get("CONDA_PREFIX")
        exe = os.environ.get("CONDA_EXE")
        conda_env = os.path.basename(os.environ.get("CONDA_DEFAULT_ENV"))
        # print(f'project:{project_dir}\tprefix:{conda_prefix}\tenv:{conda_env}\texe:{exe}')
        # cmd = 'ls'
        if conda_prefix:
            conda = f'{conda_prefix}/bin/conda'
            conda_exe = '/apps/external/anaconda/2025.12/bin/conda'
            # print(f'conda:{conda}')
            # Build the command using 'conda run -p <path_to_env>'
            cmd = [conda_exe, "run", "-n", conda_env] + cmd.split()
            # print(f'cmd:{cmd}')

        result = subprocess.Popen(cmd,
                                  env=env,
                                  stdout=self.log['stdout'], stderr=self.log['stderr'])

        return result

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
            # completed commands have already been removed by fast forward, you don't have to check
            entry = commandlist.pop()
            thiscommand = entry['command']
            self.jobid += 1
            self.log.add('main', f'Executor\tstarted\t'
                                 f'{entry["commandname"]}\t{thiscommand}\tjobid: {self.jobid}')
            job = (self.conda_run(thiscommand))
            self.joblist.append([self.jobid, job, entry['commandname'], thiscommand])
            # commandlist.remove(entry)
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
        print(f'polling: {len(self.joblist)} jobs {Log.logtime()}')
        while len(to_remove) == 0:
            for j in self.joblist:
                jid, job, stage, command = j
                # print(f'\tjob {jid} ...', end='')
                result = job.poll()
                if result is None:
                    # None indicates job is still running
                    # print(f'job {jid} still running')
                    pass

                else:
                    # job finished
                    self.running -= 1
                    self.finished += 1
                    if result == 0:
                        # success
                        self.log.add('main', f'Executor\tfinished\t{stage}\t{command}\tjobid:{jid}')
                        self.succeeded += 1

                    else:
                        # error
                        # self.log.add('stderr', f'Executor: fail jobid:{jid} stage:{stage}, ')
                        self.log.add('main',
                                     f'Executor\tjobid:{jid} failed with status={result}\t{stage}\t'
                                     f'{command}')
                        self.failed += 1

                    # include the result in the remove list, it can't be removed here because it
                    # changes the list (self.joblist) that is iterating
                    to_remove.append(j)

        # remove all finished jobs. Couldn't do it above because it shortens the joblist and
        # some jobs don't get polled
        for j in to_remove:
            self.joblist.remove(j)

        self.log.add('main', f'Executor\trunning:{self.running} succeeded:{self.finished} failed:{self.failed}')

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

    @staticmethod
    def logtime():
        """-----------------------------------------------------------------------------------------
        Create a time string for use in logs year month day hour min sec concatenated

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
        log_message = f'{Log.logtime()}\t{message}\n'
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
    w.prepare_project()

    w.command.generate()
    exec = Executor(w.command, w.log, jobs=w.option["jobs"], delay=5)
    while exec.commandlist.commands or exec.joblist:
        # continue as long as there are commands in the command list
        exec.prioritize()
        exec.startjobs()
        exec.polljobs()

        # generate more commands if possible
        w.command.generate()

    w.log.add('main', f'Project: finished {w.project}')

    exit(0)
