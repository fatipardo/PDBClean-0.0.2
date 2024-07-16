#!/usr/bin/env python
# coding: utf-8
#
from __future__ import print_function
from __future__ import division
import sys, glob
from PDBClean import pdbcleanchainstandardizationutils as chainstd

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
Structure_Sequences, structid_list, chid_list = chainstd.pdb_to_structurelists(filelist)
Standard_Sequences = {}


############################################
# INTERACTIVE ChainID STANDARDIZATION MENU #
############################################
input_menu = ""
input_menu_check_1 = ""

while(input_menu != "QUIT"):
    print("PDBClean ChainID Standardization Menu",
          "    Select one of the following options to proceed:",
          "    1) Select Standard Sequences from a chosen input structure",
          "    2) Generate Standard Sequences based on all the input structures",
          sep="\n")
    if(input_menu_check_1 == "1"):
        print("    3) Inspect/Edit Standard Sequences",
              "    4) Perform Standardization of Chain IDs",
              sep="\n")
    input_menu = input('Option Number: ')
    if (input_menu == "1"):
        Standard_Sequences, input_menu_check_1 = chainstd.select_standard_seq_from_reference(Structure_Sequences,
                                                                                             Standard_Sequences,
                                                                                             structid_list,
                                                                                             input_menu_check_1)
    elif (input_menu == "2"):
        Standard_Sequences, input_menu_check_1 = chainstd.create_standard_seq_from_consensus(Structure_Sequences,
                                                                                             Standard_Sequences,
                                                                                             chid_list,
                                                                                             input_menu_check_1)
        print("These are the standard sequences:")
        print(Standard_Sequences)
    elif (input_menu == "3" and input_menu_check_1 == "1"):
        chainstd.review_standard_seq(Structure_Sequences, Standard_Sequences)

    elif (input_menu == "4" and input_menu_check_1=="1"):
        chainstd.align_to_std_seq_and_save_to_disk(Structure_Sequences,
                                                   Standard_Sequences,
                                                   structid_list,
                                                   filelist,
                                                   target_dir=target_dir)
        print("Done!")
        input_menu = "QUIT"
