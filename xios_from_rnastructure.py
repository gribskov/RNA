"""=================================================================================================
Pipeline for generating XIOS structures using the program Fold in the RNAStructure package. By
specifying a range of delta deltaG, a series of results can be generated at different deltaG cutoffs
using a single CT file. a new CT file must be generated for each different window.

Steps:
    1) fold program (RNAstructure) - calculate suboptimal structures --> ct file
    run in --MFE mode and get the MFE
    2) convert the desired delta deltaG to percent using the MFE and run FOLD a second time to get
    the desired suboptimal structures
    3) use topology/CTread to merge structures in desired delta deltaG range and convert to XIOS

sample command:
xios_from_rnastructure.py -i ../data/fasta/ -f rnasep_a1.Buchnera_APS.fa -d 0.5,4.0,0.5
run xios_from_rnastructure.py -h for help with command line arguments
================================================================================================="""
import argparse
import glob
import os
import subprocess
import sys
import time
from datetime import datetime

from topology import RNAstructure


def formatter(prog):
    """---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Set up formatting for help
    :param prog:
    :return: argparse formatter class
    ---------------------------------------------------------------------------------------------"""
    return argparse.HelpFormatter(prog, max_help_position=50, width=120)


def options():
    """---------------------------------------------------------------------------------------------
    Command line options parsed with argparse

    :return: argparse Namespace (similar to a dictionary)
    ---------------------------------------------------------------------------------------------"""
    # print(sys.argv)

    commandline = argparse.ArgumentParser(
        description='Create XIOS from FastaA sequence using the RNAstructure fold program',
        formatter_class=formatter)
    commandline.add_argument('-q', '--quiet',
                             help='Run with minimal output to terminal',
                             action='store_true')
    commandline.add_argument('-i', '--indir',
                             help='Directory containing FastA files (%(default)s)',
                             default='./')
    commandline.add_argument('-c', '--ctdir',
                             help='Directory for output CT files (%(default)s)',
                             default='./ctfiles')
    commandline.add_argument('-x', '--xiosdir',
                             help='Directory for output XIOS files (%(default)s)',
                             default='./xiosfiles')
    commandline.add_argument('-f', '--fasta',
                             help='FastA filename, can be a glob (%(default)s)',
                             default='*.fa')
    commandline.add_argument('-p', '--percent',
                             help='Fold: percent cutoff for suboptimal structures (%(default)s)',
                             default=10)
    commandline.add_argument('-m', '--maximum',
                             help='Fold: maximum number of suboptimal structures (%(default)s)',
                             default=100)
    commandline.add_argument('-w', '--window',
                             help='Fold: window for suboptimal structures (%(default)s) or comma separated (%(default)s)',
                             default='4')
    commandline.add_argument('-d', '--ddG',
                             help='Mergestems: delta deltaG limit for supoptimal structures (%(default)s)'
                                  'or comma separated',
                             default='5.0')
    commandline.add_argument('-y', '--python',
                             help='Directory path for python executables such as ct2xios (%(default)s)',
                             default='/scratch/bell/mgribsko/rna/RNA/')
    commandline.add_argument('-r', '--rnastructure',
                             help='Directory path for RNAstructure executables such as Fold (%(default)s)',
                             default='/scratch/bell/mgribsko/rna/RNAstructure/')

    # mergecase currently unused
    # commandline.add_argument('-c', '--mergecase',
    #                          help='Mergestems: merging cases (rules) for combining suboptimal '
    #                               'structures (%(default)s). Not implemented for now.',
    #                          default='12345')

    # process command line arguments
    args = commandline.parse_args()

    # make sure paths end in /
    # d = args.__dict__
    for path in ('indir', 'ctdir', 'xiosdir', 'python', 'rnastructure'):
        if not args.__dict__[path].endswith('/'):
            args.__dict__[path] += '/'

    # ddG and window may be comma separated ranges
    # in the case of ddG, we expect a triple of begin, end, increment, if not specified, increment is 1.0
    # first the defaults
    if args.ddG.find(','):
        minmax = args.ddG.split(',')
        args.ddG_min = float(minmax[0])
        args.ddG_max = args.ddG_min
        args.ddG_inc = 1.0
        if len(minmax) > 1:
            # ddg_max is present
            args.ddG_max = float(minmax[1])
        if len(minmax) > 2:
            # ddg_inc is present
            args.ddG_inc = float(minmax[2])

    else:
        args.ddG_max = float(args.ddG)
        args.ddG_min = float(args.ddG)
        args.ddG_inc = 1.0

    if args.window.find(','):
        minmax = args.window.split(',')
        if len(minmax) > 1:
            [args.window_min, args.window_max] = minmax
            args.window_min = int(args.window_min)
            args.window_max = int(args.window_max)
        else:
            args.window_max = args.window_min = int(args.window)

    return args


def safe_mkdir(args, dirname):
    """---------------------------------------------------------------------------------------------
    Try to create a directory, if directory already exists, use existing directory

    :param args: namespace  command line arguments from argparse
    :param dirname: string, path to directory
    :return: True
    ---------------------------------------------------------------------------------------------"""
    try:
        os.mkdir(dirname)
        return True
    except FileExistsError:
        if not args.quiet:
            sys.stderr.write(f'Directory {dirname} already exists, using existing directory\n')
    #       exit(2)

    return True


def safe_file(filename, mode):
    """---------------------------------------------------------------------------------------------
    check if the named file is readable, mode = 'r', or does not already exist, mode=w.
    Exit with status = 3 or status = 4, respectively if the tests fail

    :param filename: string
    :param mode: string 'r' or 'w'
    :return: True
    ---------------------------------------------------------------------------------------------"""
    if mode == 'r':
        # print(f'safe_file r-branch file: {filename}')
        if not os.access(filename, os.R_OK):
            sys.stderr.write(f'File {filename} cannot be read\n')
            exit(3)

    else:
        if os.path.exists(filename):
            sys.stderr.write(f'File {filename} already exists, move or delete {filename}\n')
            return False
            # exit(4)
        # print(f"apparently {filename} doesn't exist")

    return True


def ct_from_fasta(args, fasta):
    """---------------------------------------------------------------------------------------------
    create a name for a ct file from the fastafile name.
    1. remove any directory path
    2. remove the last suffix if it is .fa, or .fasta
    3. add the suffix .ct

    :param args: Namespace, command line arguments
    :param fasta: string, name of fasta file
    :return ct: string, name of ct file
    ---------------------------------------------------------------------------------------------"""
    ct = os.path.basename(fasta)
    l = len(ct)
    if ct.rindex('.fa') == l - 3:
        ct = ct[:-3]
    elif ct.rindex('.fasta') == l - 6:
        ct = ct[:-6]

    ct += f'.w{args.window}'
    ct += '.ct'
    return f'{args.ctdir}{ct}'


def xios_from_ct(args, ct):
    """---------------------------------------------------------------------------------------------
    create a name for a xios file from the ct name.
    1. remove any directory path
    2. remove the last suffix if it is .ct
    3. add the suffix .xios

    :param args: Namespace, command line arguments, contains ddG
    :param ct: string, name of ct file
    :return: string, name of xios file
    ---------------------------------------------------------------------------------------------"""
    xios = os.path.basename(ct)
    l = len(xios)
    if xios.rindex('.ct') == l - 3:
        xios = xios[:-3]
    xios += f'.d{args.ddg}'
    # xios += f'.c{args.mergecase}'

    xios += '.xios'

    return xios


def get_mfe_from_ct(CT):
    """---------------------------------------------------------------------------------------------
    Read the first energy from the CT file the first energy is the MFE
    396  ENERGY = -27.6  tmRNA.Gracilaria_tenuistipitata.AY673996.fa

    :param CT: string, filename of input CT file
    :return: float, MFE
    ---------------------------------------------------------------------------------------------"""
    ct = None
    if safe_file(CT, 'r'):
        ct = open(CT, 'r')

    field = []
    for line in ct:
        if line.find('ENERGY'):
            field = line.split()
            break

    ct.close()
    return float(field[3])


def runfold(args, fasta, ct, percent):
    """---------------------------------------------------------------------------------------------
    Fold is run twice, the first time to just get the MFE. The second time under the ddG and window
    parameters to get the suboptimal alignments

    Run the RNAstructure Fold program
    USAGE: Fold <sequence file> <CT file> [options]
    All flags are case-insensitive, and grouping of flags is not allowed.

    ==== Required Parameters ====
    <sequence file>
        The name of a file containing an input sequence. Acceptable formats include
        SEQ, FASTA and raw-sequence plain-text files.
        If the name is a hyphen (-), the file will be read from standard input
        (STDIN).

    <CT file>
        The name of a CT file to which output will be written. If the --bracket
        flag is present, output will be written as a dot-bracket file.
        If the file name is a hyphen (-), the structure will be written to standard
        output (STDOUT) instead of a file.

    ==== Selected Options ====

    -l -L --loop
        Specify a maximum internal/bulge loop size.
        Default is 30 unpaired numcleotides.

    -m -M --maximum
        Specify a maximum number of structures.
        Default is 20 structures.

    --name
        Specify a name for the sequence. This will override the name in the
        sequence file.

    -p -P --percent
        Specify a maximum percent energy difference.
        Default is 10 percent (specified as 10, not 0.1).

    -s -S --save
        Specify the name of a save file, needed for dotplots and refolding.
        Default is not to generate a save file.

    -t -T --temperature
        Specify the temperature at which calculation takes place in Kelvin.
        Default is 310.15 K, which is 37 degrees C.
    -w -W --window
        Specify a window size.
        Default is determined by the length of the sequence.

    --warnings --warn
        Set the behavior for non-critical warnings (e.g. related to invalid
        nucleotide positions or duplicate data points in SHAPE data). Valid values
        are:
          * on  -- Warnings are written to standard output. (default)
          * err -- Warnings are sent to STDERR. This can be used in automated scripts
            etc to detect problems.
          * off -- Do not display warnings at all (not recommended).

    :param args: Namespace  command line arguments from options() using argparse
    :param fasta: string    readable fasta file
    :param ct: string       writable ct file
    :return: string         comment for inclusion in XIOS file
    ---------------------------------------------------------------------------------------------"""

    # check the input (fa) and output files. Fail if input file is not available/readable.
    # if ct file is present, assume this sequence has been processed already and skip
    safe_file(fasta, 'r')
    if safe_file(ct, 'w'):
        # ct file does not exist
        # print(f'ctfile {ct} is writable')
        pass
    else:
        # probably ctfile exists so skip
        # print(f'ct file {ct} exists. {fasta} skipped')
        return f'{fasta} skipped'
    # -----------------------------------------------------------------------------------------------
    # pass 1 to get DeltaG
    # -----------------------------------------------------------------------------------------------
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    comment = f'\tRNAstructure/Fold {now}\n'

    if not args.quiet:
        print(f'\trunning fold pass 1: {fasta} => {ct}')
    exe = args.rnastructure + '/exe/Fold'
    opt = [exe, fasta, ct]
    opt += ['-mfe']
    opt += ['-t', '298']
    opt += ['-l', '50', '-y']
    result = subprocess.run(opt, capture_output=True)

    mfe = get_mfe_from_ct(ct)
    comment += f'\tMFE: {mfe} delta(deltaG): {args.ddG_max}\n'
    try:
        percent = int(-100 * args.ddG_max / mfe) + 1  # round up to make sure you get all structures
    except ZeroDivisionError:
        percent = args.percent

    # -----------------------------------------------------------------------------------------------
    # pass 2 to get suboptimal folds
    # -----------------------------------------------------------------------------------------------
    if not args.quiet:
        print(f'\trunning fold pass 2: {fasta} => {ct}')
    exe = args.rnastructure + '/exe/Fold'
    opt = [exe, fasta, ct]
    opt += ['-p', f'{percent}']
    opt += ['-m', f'{args.maximum}']
    opt += ['-w', f'{args.window}']
    opt += ['-t', '298']
    opt += ['-l', '50', '-y']
    result = subprocess.run(opt, capture_output=True)

    comment += f'\tcommand options: {" ".join(opt[3:])}\n'
    comment += f'\tinput: {fasta}\n'
    comment += f'\toutput: {ct}\n'

    return comment


def runmergestems(arg, ct, xios, comment):
    """---------------------------------------------------------------------------------------------
    CT files are converted using CTread in Topology.
    TODO since the only thing passed in args is the ddG, it'd be more transparent to use an explicit
    argument
    :param arg: namespace   command line argument namespace, contains ddG_max
    :param ct: string       readable ct file
    :param xios: string     name of xios file to generate
    :param comment: string  string to add to XIOS comment block when topology is created in
                            runmergestructure()
    :return: string         updated comment string for later use
    ---------------------------------------------------------------------------------------------"""
    if safe_file(xios, 'w'):
        pass
    else:
        return f'runmergestems skipping {xios}'

    xiosout = open(xios, 'w')

    rna = RNAstructure()
    # add the comment passed in from runfold
    rna.comment.append(comment)
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    comment = f'\trunmergestems delta(deltaG) = {args.ddg} {now}\n'
    comment += f'\tinput_file: {ct}\n'
    comment += f'\toutput_file: {xios}\n'
    rna.comment.append(comment)

    rna.CTRead(ct, ddG=args.ddg)
    rna.adjacency_from_stemlist()
    rna.edgelist_from_adjacency(include="ijo", whole=False)

    rna.XIOSwrite(xiosout)

    return


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    now = time.localtime()
    args = options()

    # code locations, see options() to alter defaults
    os.environ['DATAPATH'] = args.rnastructure + 'data_tables'
    input = args.indir + args.fasta

    if args.quiet:
        print(f'xios_from_rnastructure {time.asctime(now)}', end='')
        print(f'\tdirs;fasta:{args.indir};ct:{args.ctdir};xios:{args.xiosdir}\tinput:{input}')
    else:
        print('Paths')
        print(f'\tRNAstructure DATAPATH = {os.environ["DATAPATH"]}')
        print(f'\tRNAstructure executables: {args.rnastructure}')
        print(f'\tPython executables: {args.python}')

        print('Parameters')
        print(f'\tFold:percent={args.percent}')
        print(f'\tFold:maximum={args.maximum}')
        print(f'\tFold:window={args.window_min},{args.window_max}')
        print(f'\tdelta deltaG={args.ddG_min},{args.ddG_max},{args.ddG_inc}')
        # print(f'\tmergstems:mergecases={args.mergecase}')

        print('Inputs and outputs')
        print(f'\tFasta files: {args.indir}')
        print(f'\tCT files: {args.ctdir}')
        print(f'\tXIOS files: {args.xiosdir}\n')

    # check if there are inputs and create directories

    fastafiles = glob.glob(input)
    if not fastafiles:
        sys.stderr.write('No files match the specified input ({})\n'.format(input))
        exit(1)

    safe_mkdir(args, args.ctdir)
    safe_mkdir(args, args.xiosdir)

    # for each FastA file, generate the CT file, then convert to XIOS
    commentfold = ''
    ctlist = []
    fa_n = 0
    ct_n = 0
    xios_n = 0
    for fasta in fastafiles:
        fa_n += 1
        if not args.quiet:
            print(f'\nprocessing {fasta}')

        for window in range(args.window_min, args.window_max + 1):
            # run fold for each window size, CT files go to args.ctdir
            args.window = window
            ct = ct_from_fasta(args, fasta)
            commentfold = runfold(args, fasta, ct, percent=0)

            # all XIOS output goes to args.xiosdir
            # the stems that are the same in structures with different ddG are merged

            ddG = args.ddG_min
            while ddG <= args.ddG_max:
                # print(f'ddg loop')
                args.ddg = ddG
                xios = f'{args.xiosdir}/{xios_from_ct(args, ct)}'
                # sys.stderr.write(f'filecheck ct={ct}\txios={xios}\n')
                runmergestems(args, ct, xios, commentfold)
                xios_n += 1
                ddG += args.ddG_inc

    # final report
    if not args.quiet:
        print(f'\nFastA files processed: {fa_n}')
        print(f'CT files processed: {ct_n}')
        print(f'XIOS files produced: {xios_n}')

    exit(0)
