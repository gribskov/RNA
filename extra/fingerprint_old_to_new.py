"""=====================================================================================================================
fingerprint_old_to_new.py
Convert old XML format to current YAML format. should be lossless

usage:
    python fingerprint_old_to_new.py in.xpt
    output is the same as the input with xpt changed to fpt

2026-03-30 gribskov
====================================================================================================================="""
import sys
import os
import xios
import topology
from fingerprint import Fingerprint

# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    xptfile = sys.argv[1]
    xpt = Fingerprint()
    xpt.readXML(xptfile)

    xpt.information['fingerprint']['converted from'] = f'{xptfile.rstrip()} at: {xpt.setdate()}'
    fptfile = os.path.basename(xptfile).replace('.xpt', '.fpt')
    xpt.writeYAML(fptfile, sort='alpha')

    exit(0)
