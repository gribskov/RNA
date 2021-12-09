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


def formatter(prog):
    """---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Set up formatting for help
    :param prog:
    :return: argparse formatter class
    ---------------------------------------------------------------------------------------------"""
    return argparse.HelpFormatter(prog, max_help_position=50, width=100)


def options():
    """---------------------------------------------------------------------------------------------
    Command line options parsed with argparse

    :return: argparse Namespace (similar to a dictionary)
    ---------------------------------------------------------------------------------------------"""
    commandline = argparse.ArgumentParser(
        description='Create XIOS from FastaA sequence using the RNAstructure fold program',
        formatter_class=formatter)
    commandline.add_argument('-i', '--indir',
                             help='Directory containing FastA files',
                             default='./')
    commandline.add_argument('-f', '--fasta',
                             help='FastA filename, can be a glob (%(default)s)',
                             default='*.fa')
    commandline.add_argument('-p', '--percent',
                             help='Fold: percent cutoff for suboptimal structures (%(default)s)',
                             default=4)
    commandline.add_argument('-m', '--maximum',
                             help='Fold: maximum number of suboptimal structures (%(default)s)',
                             default=20)
    commandline.add_argument('-w', '--window',
                             help='Fold: window for suboptimal structures (%(default)s)',
                             default=4)
    commandline.add_argument('-d', '--ddG',
                             help='Mergestems: delta deltaG limit for supoptimal structures (%(default)s)',
                             default=5)
    commandline.add_argument('-c', '--mergecase',
                             help='Mergestems: merging cases (rules) for combining suboptimal '
                                  'structures (%(default)s)',
                             default='12345')

    # process command line arguments
    args = commandline.parse_args()
    if not args.indir.endswith('/'):
        args.indir += '/'

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
        sys.stderr.write(f'Directory {dirname} already exists, move or delete {dirname}')
        exit(2)


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


def ct_from_fasta(fasta):
    """---------------------------------------------------------------------------------------------
    create a name for a ct file from the fastafile name.
    1. remove any directory path
    2. remove the last suffix if it is .fa, or .fasta
    3. add the suffix .ct

    :param fasta: string, name of fasta file
    :return ct: string, name of ct file
    ---------------------------------------------------------------------------------------------"""
    ct = os.path.basename(fasta)
    l = len(ct)
    if ct.rindex('.fa') == l - 3:
        ct = ct[:-3]
    elif ct.rindex('.fasta') == l - 6:
        ct = ct[:-6]

    ct += '.ct'
    return ct


def xios_from_ct(ct):
    """---------------------------------------------------------------------------------------------
    create a name for a xios file from the ct name.
    1. remove any directory path
    2. remove the last suffix if it is .ct
    3. add the suffix .xios

    :param ct: string, name of ct file
    :return ct: string, name of xios file
    ---------------------------------------------------------------------------------------------"""
    xios = os.path.basename(ct)
    l = len(xios)
    if xios.rindex('.ct') == l - 3:
        xios = xios[:-3]

    xios += '.xios'
    return xios


def runfold(args, fasta, ct):
    """---------------------------------------------------------------------------------------------
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
    :return:
    ---------------------------------------------------------------------------------------------"""
    safe_file(fasta, 'r')
    safe_file(ct, 'w')

    exe = args.rnastructure + '/Fold'
    opt = [exe, fasta, ct]
    opt += ['-p', f'{args.percent}']
    opt += ['-m', f'{args.maximum}']
    subprocess.call(opt)

    return


def runmergestems(mergecases, ddG, ct):
    """---------------------------------------------------------------------------------------------
    mergestem.pl [-c <mergecases>] [-g <ddG>] [l <int>] [-p <plot_file>]  [-q] [-m <curated_file>] 
            <CT_file>
    
    Mergestem compares that multiple structures in the CT file and convertss them to a single XIOS 
    structure.
    
    :param mergecases: string, list of mergecases to apply
    :param ddG: maximum delta delta G for suboptimal structures
    :param ct: string, readable ct file
    :return: 
    ---------------------------------------------------------------------------------------------"""
    safe_file(xios, 'w')

    exe = args.perl_src + '/mergestems.pl'
    opt = [exe, ct]
    opt += ['-c', f'{args.mergecase}']
    opt += ['-g', f'{args.ddg}']
    subprocess.call(opt)

    return


# --------------------------------------------------------------------------------------------------
# main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    args = options()

    # code locations, add to args
    rnastructure = "/scratch/bell/mgribsko/rna/RNAstructure/exe"
    args.rnastructure = rnastructure
    perl_src = "/depot/mgribsko/rna/RNA/perl_src"
    args.perl_src = perl_src
    print(f'RNAstructure executables: {rnastructure}')
    print(f'Mergestem executable: {perl_src}')

    # check if there are inputs and create directories

    input = args.indir + args.fasta
    print(f'FastaA input: {input}')
    fastafiles = glob.glob(input)

    if not fastafiles:
        sys.stderr.write('No files match the specified input ({})\n'.format(input))
        exit(1)

    ctdir = 'ctfiles'
    safe_mkdir(ctdir)
    xiosdir = 'xiosfiles'
    safe_mkdir(xiosdir)

    for fasta in fastafiles:
        # run fold
        ct = 'ctfiles/' + ct_from_fasta(fasta)
        runfold(args, fasta, ct)

        # run mergestems
        xios = 'xiosfiles' + xios_from_ct(ct)
        runmergestems(args, ct)

    exit(0)
