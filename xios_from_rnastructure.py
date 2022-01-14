"""=================================================================================================
Pipeline for generating XIOS structures using the program fold in the RNAStructure package.

Steps:
    fold program (RNAstructure) - calculate suboptimal structures --> ct file
    mergestems.pl -- merge multiple suboptimal structures --> XIOS file

================================================================================================="""
import argparse
import glob
import os
import subprocess
import sys
import time

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
                             help='Fold: window for suboptimal structures (%(default)s) or comma separated',
                             default='4')
    commandline.add_argument('-d', '--ddG',
                             help='Mergestems: delta deltaG limit for supoptimal structures (%(default)s) '
                                  'or comma separated',
                             default='5')
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
    if args.ddG.find(','):
        minmax = args.ddG.split(',')
        if len(minmax) > 1:
            [args.ddG_min, args.ddG_max] = minmax
            args.ddG_min = int(args.ddG_min)
            args.ddG_max = int(args.ddG_max)
        else:
            args.ddG_max = args.ddG_min = int(args.ddG)

    if args.window.find(','):
        minmax = args.window.split(',')
        if len(minmax) > 1:
            [args.window_min, args.window_max] = minmax
            args.window_min = int(args.window_min)
            args.window_max = int(args.window_max)
        else:
            args.window_max = args.window_min = int(args.window)

    return args


def safe_mkdir(dirname):
    """---------------------------------------------------------------------------------------------
    Try to create a directory, if directory already exists, exit with status = 2

    :param dirname: string, path to directory
    ---------------------------------------------------------------------------------------------"""
    try:
        os.mkdir(dirname)
        return True
    except FileExistsError:
        sys.stderr.write(f'Directory {dirname} already exists, using existing directory')
    #    exit(2)

    return True


def safe_file(filename, mode):
    """---------------------------------------------------------------------------------------------
    check if the named file is readable, mode = 'r', or does not already exist, mode=w.
    Exit with status = 3 or status = 4, respectively if the tests fail

    :param filename: string
    :param mode: string, 'r' or 'w'
    :return:
    ---------------------------------------------------------------------------------------------"""
    if mode is 'r':
        if not os.access(filename, os.R_OK):
            sys.stderr.write(f'File {filename} cannot be read')
            exit(3)

    else:
        if os.path.exists(filename):
            sys.stderr.write(f'File {filename} already exists, move or delete {filename}')
            exit(4)

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
    return ct


def xios_from_ct(args, ct):
    """---------------------------------------------------------------------------------------------
    create a name for a xios file from the ct name.
    1. remove any directory path
    2. remove the last suffix if it is .ct
    3. add the suffix .xios

    :param args: Namespace, command line arguments
    :param ct: string, name of ct file
    :param ddG: int, delta deltaG cutoff value
    :return ct: string, name of xios file
    ---------------------------------------------------------------------------------------------"""
    xios = os.path.basename(ct)
    l = len(xios)
    if xios.rindex('.ct') == l - 3:
        xios = xios[:-3]
    xios += f'.d{args.ddg}'
    # xios += f'.c{args.mergecase}'

    xios += '.xios'
    print(f'ddg file:{xios}')
    return xios


def get_mfe_from_ct(CT):
    """---------------------------------------------------------------------------------------------
    Read the first energy from the CT file the first energy is the MFE

    :param ct: string, filename
    :return: float, MFE
    ---------------------------------------------------------------------------------------------"""
    ct = safe_file(CT, 'r')
    for line in ct:
        if line.find('ENERGY'):
            field = line.split()

    ct.close()

    return field[3]


def runfold(args, fasta, ct, percent):
    """---------------------------------------------------------------------------------------------
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

    :param args: Namespace, command line arguments from options()
    :param fasta: string, readable fasta file
    :param ct: string, writable ct file
    :return: string, comment for inclusion in XIOS file
    ---------------------------------------------------------------------------------------------"""
    safe_file(fasta, 'r')
    safe_file(ct, 'w')

    exe = args.rnastructure + '/exe/Fold'
    opt = [exe, fasta, ct]
    opt += ['-p', f'{percent}']
    opt += ['-m', f'{args.maximum}']
    opt += ['-w', f'{args.window}']
    subprocess.call(opt)

    comment = f'{exe} -p {args.percent} -m {args.maximum} -w {args.window} {fasta} {ct}\n'
    comment += f'{time.asctime(time.localtime())}'

    return comment


def runmergestems(arg, ct, xios, comment):
    """---------------------------------------------------------------------------------------------
    mergestem.pl [-c <mergecases>] [-g <ddG>] [l <int>] [-p <plot_file>]  [-q] [-m <curated_file>] 
            <CT_file>
    
    Mergestem compares that multiple structures in the CT file and converts them to a single XIOS
    structure.
    
    :param mergecases: string, list of mergecases to apply
    :param ddG: maximum delta delta G for suboptimal structures
    :param ct: string, readable ct file
    :param commen: string, string to add to XIOS comment block
    :return: 
    ---------------------------------------------------------------------------------------------"""
    safe_file(xios, 'w')

    xiosout = open(xios, 'w')

    exe = args.python + '/ct2xios.py'
    opt = ['python3', exe, ct]
    # opt += ['-c', f'{args.mergecase}']
    opt += ['-d', f'{args.ddg}']
    subprocess.call(opt, stdout=xiosout)
    xiosout.close()

    rna = RNAstructure()
    rna.comment.append(comment)
    rna.CTRead(ct, ddG=args.ddg)
    # now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    # rna.comment.append('creation_date {}'.format(now))
    # rna.comment.append('input_file {}'.format(ctfile))
    # if rna.energy:
    #     rna.comment.append('input_format unifold')
    # else:
    #     rna.comment.append('input_format RNAstructure')
    rna.adjacency_from_stemlist()
    rna.edgelist_from_adjacency(include="ijo", whole=False)

    rna.XIOSwrite(sys.stdout)

    return


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    now = time.localtime()
    args = options()

    # code locations, see options() to alter defaults
    os.environ['DATAPATH'] = args.rnastructure + 'data_tables'
    print(f'xios_from_rnastructure {time.asctime(now)}\n')
    if not args.quiet:
        print('Paths')
        print(f'\tRNAstructure DATAPATH = {os.environ["DATAPATH"]}')
        print(f'\tRNAstructure executables: {args.rnastructure}')
        print(f'\tPython executables: {args.python}')

    print('Parameters')
    print(f'\tFold:percent={args.percent}')
    print(f'\tFold:maximum={args.maximum}')
    print(f'\tFold:window={args.window}')
    print(f'\tdelta deltaG={args.ddG}')
    # print(f'\tmergstems:mergecases={args.mergecase}')

    input = args.indir + args.fasta
    print('Inputs and outputs')
    print(f'\tFasta files: {input}')
    print(f'\tCT files: {args.ctdir}')
    print(f'\tXIOS files: {args.xiosdir}')

    # check if there are inputs and create directories

    fastafiles = glob.glob(input)
    if not fastafiles:
        sys.stderr.write('No files match the specified input ({})\n'.format(input))
        exit(1)

    safe_mkdir(args.ctdir)
    safe_mkdir(args.xiosdir)

    # for each FastA file, generate the CT file, then convert to XIOS
    commentfold = {}
    ctlist = []
    fa_n = 0
    ct_n = 0
    xios_n = 0
    for fasta in fastafiles:
        fa_n += 1
        if not args.quiet:
            print(f'processing {fasta}')

        for window in range(args.window_min, args.window_max + 1):
            # run fold for each window size, CT files go to args.ctdir
            args.window = window
            ct = 'ctfiles/' + ct_from_fasta(args, fasta)
            ctlist.append(ct)
            commentfold[ct] = runfold(args, fasta, ct, percent=0)

            mfe = get_mfe_from_ct(ct)
            try:
                percent = int(100 * args.ddG / mfe)
            except ZeroDivisionError:
                percent = args.percent

            commentfold[ct] = runfold(args, fasta, ct, percent=percent)
            ct_n += 1

            # all XIOS output goes to the args.xiosdir
            # for each ct file we can use different ddG
            # the stems that are the same in structures with different ddG are merged

            for ddG in range(args.ddG_min, args.ddG_max + 1):
                args.ddg = int(ddG)
                xios = f'{args.xiosdir}/{xios_from_ct(args, ct)}'
                # sys.stderr.write(f'filecheck ct={ct}\txios={xios}\n')
                runmergestems(args, ct, xios, commentfold[ct])
                xios_n += 1

    # final report
    print(f'FastA files processed: {fa_n}')
    print(f'CT files processed: {ct_n}')
    print(f'XIOS files produced: {xios_n}')

    exit(0)
