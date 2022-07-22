#!/Users/fatima/anaconda3/envs/PDBC0722/bin/python
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
input_menu_check_2 = ""

while(input_menu != "QUIT"):
    print("PDBClean ChainID Standardization Menu",
          "    Select one of the following options to proceed:",
          "    1) Select Standard Sequences from input structure",
          "    2) Create Standard Sequences from consensus of input structures",
          sep="\n")
    if(input_menu_check_1 == "1"):
        print("    3) Inspect/Edit Standard Sequences",
              "    4) Perform pairwise alignments against Standard Sequences",
              sep="\n")
    if(input_menu_check_2 == "1"):
        print("    5) Inspect/Edit chain ID reassignments",
              "    6) Perform Standardization of Chain IDs",
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
    elif (input_menu == "3" and input_menu_check_1 == "1"):
        chainstd.review_standard_seq(Structure_Sequences, Standard_Sequences)
    elif (input_menu == "4" and input_menu_check_1 == "1"):
        ChainReassignmentMapping_List, ChainReassignmentScores_List, input_menu_check_2 = chainstd.align_to_standard_seq(Structure_Sequences, 
                                                                                                                         Standard_Sequences,
                                                                                                                         structid_list)
    elif (input_menu == "6" and input_menu_check_2 == "1"):
        chainstd.reassignedmaps_to_pdb(filelist, ChainReassignmentMapping_List, structid_list, target_dir=target_dir)
        chainstd.reassignedmaps_to_log(ChainReassignmentMapping_List, ChainReassignmentScores_List, structid_list, target_dir=target_dir)
        print("Done!")
        input_menu="QUIT"
