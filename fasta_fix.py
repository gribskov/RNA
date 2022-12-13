"""=================================================================================================
Check fasta files and make sure that all bases are uppercase and ACGT

Usage:
    fasta_fix.py <input_directory> <output_directory>
================================================================================================="""
import os
import sys


def read_directories(path, dir_list):
    """=================================================================================================================
    This function accepts a path and list of directories inside path. It reads each file in each directory and writes
    into a file any unacceptable characters.
    :param path: full directory path
    :param dir_list: list of subdirectories
    ================================================================================================================="""

    bases = ['A', 'C', 'T', 'G']
    fa_dict = dict()
    for dir in dir_list:  # for each directory in dir_list
        dir_path = path + dir
        files = os.listdir(dir_path)
        for file in files:  # for each file in subdirectory
            split_file = file.split('.')
            if split_file[1] == 'fa' or split_file[1] == 'fasta':  # make sure file is fasta file
                file_path = (dir_path + '/' + file)
                with open(file_path, 'r') as fa_file:  # open fasta file
                    for line in fa_file:  # for each line in fasta file
                        line = line.rstrip()
                        if line.startswith('>'):  # skip lines that start with '>'
                            continue
                        for letter in line:  # for each letter in line
                            if letter not in bases:  # check letter against those in bases list
                                if file_path not in fa_dict.keys():  # add file path to dictionary keys
                                    fa_dict.update({file_path: str()})
                                fa_dict[
                                    file_path] += letter  # add letter to dictionary using file path as key
                                print(f'{file_path} had an oopsie')

    # Write out dictionary keys and values
    with open('bad_bases.txt', 'w') as file:
        for key in fa_dict.keys():
            file.write(key + '\n')
            file.write(fa_dict[key] + '\n')


def fasta_read(filename):
    """---------------------------------------------------------------------------------------------
    read the fasta file and return as a dictionary = {'id', 'documentation', 'sequence'}

    :param filename: string, readable input file
    :return: dict, sequence dictionary
    ---------------------------------------------------------------------------------------------"""
    try:
        fasta_in = open(filename, 'r')
    except OSError:
        sys.stderr.write(f'fasta_fix:fasta_read - cannot open fasta file ({filename})')
        exit(1)

    fasta = {'id':'', 'documentation':'', 'sequence':''}
    for line in fasta_in:
        if line.lstrip().startswith('>'):
            doc = ''
            id = ''
            try:
                id, doc = line.strip().split()
            except ValueError:
                # no documentation
                id = line.strip()

            fasta['id'] = id.replace('>', '')
            fasta['documentation'] = doc
            fasta['sequence'] = ''
        else:
            fasta['sequence'] += line.strip()

    fasta_in.close()

    return fasta


def fasta_write(filename, fasta, linelen=100):
    """---------------------------------------------------------------------------------------------
    write the sequence dictionary = {'id', 'documentation', 'sequence'}

    :param filename: string, openable output file
    :param fasta: dict, {'id', 'documentation', 'sequence'}
    :return: True
    ---------------------------------------------------------------------------------------------"""
    try:
        fasta_out = open(filename, 'w')
    except OSError:
        sys.stderr.write(f'fasta_fix:fasta_write - cannot open fasta file ({filename})')
        exit(2)

    fasta_out.write(f">{fasta['id']} {fasta['documentation']}\n")
    begin = 0
    while begin < len(fasta['sequence']):
        fasta_out.write(f"{fasta['sequence'][begin:begin + linelen]}\n")
        begin += linelen

    fasta_out.close()
    return True


def fasta_fixbases(fasta):
    """---------------------------------------------------------------------------------------------
    Sequences may contain not ACGT bases, convert them to A

    :paraam fasta: dict, , {'id', 'documentation', 'sequence'}
    :return: fasta, int, number of bases changed
    ---------------------------------------------------------------------------------------------"""
    newbase = 'A'
    base = list(fasta['sequence'])
    changed = 0
    for i in range(len(base)):
        if base[i] not in 'ACGTU':
            base[i] = newbase
            changed += 1

    return ''.join(base), changed


# ==================================================================================================
# main/test
# ==================================================================================================
if __name__ == '__main__':

    fasta_suffix = '.fa'

    indir = sys.argv[1]
    if not indir.endswith('/'):
        indir += '/'
    outdir = sys.argv[2]
    if not outdir.endswith('/'):
        outdir += '/'

    for file in os.listdir(indir):
        # for each file input directory
        if file.endswith(fasta_suffix):
            fa = fasta_read(indir + file)
            fa_original = fa['sequence'][:]
            fa['sequence'] = fa['sequence'].upper()
            if fa['sequence'] != fa_original:
                fa['documentation'] += " original bases converted to uppercase;"
            fixed, changed = fasta_fixbases(fa)
            if changed:
                fa['documentation'] += f' {changed} original bases changed to A;'
                fa['sequence'] = fixed
            fasta_write(outdir + file, fa)

    exit(0)
