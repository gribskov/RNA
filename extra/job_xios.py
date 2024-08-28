#
# Ben Iovino
# BIOL494, 1/21/22
# Job
# This program looks through every fasta file in a directory and uses it
# as an input for a program in the RNAstructure package.
#


# Import libraries
import os
import subprocess


# Define a function that goes through a directory and makes a list of each
# file in the directory.
def get_file_list(directory):
    # print(f'file dir = {directory}')
    # Initialize a list of file names
    filelist = list()

    # Read each file name in the directory path
    for file in os.listdir(directory):
        # print(f'file={file}')
        # Append file name to list
        # filelist.append(directory + "/" + file)
        filelist.append(file)

    # Return list of file names
    return filelist


# Define a function that performs a job on each file in the filelist
def do_job(filelist, directory):
    # Use a for loop to iterate through each file in filelist
    for file in filelist:
        # Use subprocess to run Unix command
        # exe = 'python /scratch/scholar/biovino/xios_from_rnastructure.py'
        command = ['python', '/scratch/bell/mgribsko/rna/RNA/xios_from_rnastructure.py']
        command += ['-i', directory]
        command += ['-f', file]
        command += ['-c', './ctdir']
        command += ['-x', './xiosdir']
        command += ['-w', '4']
        command += ['-d', '5']

        subprocess.call(command)


# Define the main function to ask for a directory path and run the job function
def main():
    # Enter directory path
    # directory = "/scratch/scholar/biovino/RNAdata/Test"
    directory = '/scratch/bell/mgribsko/rna/testjob'
    # print(f'dir={directory}')

    # Call the get_file_list function providing the directory path as a parameter
    filelist = get_file_list(directory)
    # print(f'filelist: {filelist}')

    # Call do_job function
    do_job(filelist, directory)


main()
exit(0)
# End the main function
