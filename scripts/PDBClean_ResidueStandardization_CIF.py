#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
from __future__ import division
import sys, glob
from PDBClean import pdbcleanresiduestandardizationutils as resstd

########################
# READ INPUT ARGUMENTS #
########################
n_arg = len(sys.argv)
if(n_arg<3):
    print('Usage error: {0} <source directory> <target directory>'.format(sys.argv[0]))
    sys.exit()
source_dir=sys.argv[1]
target_dir=sys.argv[2]


#############################################
# READ PDB FILES AND DEFINE STRUCTURE LISTS #
#############################################
filelist=glob.glob(source_dir+'/*.cif')
Structure_Sequences, ChID_ResiNum_Vector, structid_list, chid_list = resstd.pdb_to_structurelists(filelist)


############################################
# INTERACTIVE RESIDUE STANDARDIZATION MENU #
############################################
input_menu = ""
input_menu_check = ""

while(input_menu != "QUIT"):
    print("PDBClean Residue Number Standardization Menu",
          "    Select one of the following options to proceed:",
          "    1) Perform multiple alignments to identify residues",
          sep="\n")
    if(input_menu_check == "1"):
        print("    2) View conversion template")
    if(input_menu_check == "1"):
        print("    3) Perform residue number standardization")
    input_menu = input('Option Number: ')
    if (input_menu == "1"):
        Structure_Sequences_Aligned, Structure_ConversionTemplate, chid_list, input_menu_check = resstd.perform_multiple_alignment(Structure_Sequences, 
                                                                                                                                   ChID_ResiNum_Vector, 
                                                                                                                                   structid_list, 
                                                                                                                                   chid_list, 
                                                                                                                                   input_menu_check)
    elif (input_menu == "2" and input_menu_check == "1"):
        resstd.show_conversiontemplate(Structure_ConversionTemplate)
    elif (input_menu == "3" and input_menu_check == "1"):
        resstd.conversiontemplate_to_pdb(filelist, Structure_ConversionTemplate, target_dir=target_dir)
