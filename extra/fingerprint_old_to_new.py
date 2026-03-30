"""=====================================================================================================================
fingerprint_old_to_new.py
Convert old XML format to current YAML format. should be lossless

2026-03-30 gribskov
====================================================================================================================="""
import xios
import topology
from fingerprint import Fingerprint
# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':

    xptfile = '../data/curated_fingerprint/rnasep_b.Mycoplasma_fermentans.xios.xpt'
    xpt = Fingerprint()
    xpt.readXML(xptfile)

    xpt.information['fingerprint']['converted from'] = f'{xptfile.rstrip()} at: {xpt.setdate()}'
    fptfile = xptfile.replace('.xpt', '.fpt')
    xpt.writeYAML(fptfile)

    exit(0)
