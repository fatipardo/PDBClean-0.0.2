#!/usr/bin/env python
# coding: utf-8
#
# ! ! ! master_molID_class_list is very important
# This is the list that contains every file's MolID class
# !! molIDConversion_list is important . . . it contains the objects
# MolIDConversion and is what is going to be updated by the user and evaulated
# to determine when the next step in the program is unlocked
# Create list of MolIDConversion objects using unique_molID_occur_map

from __future__ import print_function
import sys, glob
from PDBClean import pdbcleanmolidcifutils as molidutils


########################
# READ INPUT ARGUMENTS #
########################
n_arg = len(sys.argv)
if(n_arg<3):
    print('Usage error: {0} <source directory> <target directory>'.format(sys.argv[0]))
    sys.exit()
source_dir=sys.argv[1]
target_dir=sys.argv[2]


#########################################
# READ PDB FILES AND DEFINE MolID LISTS #
#########################################

filelist=glob.glob(source_dir+'/*.cif')
master_molID_class_list = molidutils.pdb_to_masterlist(filelist)
unique_molID_occur_map  = molidutils.CreateMasterUniqueMolIDMap(master_molID_class_list)
molIDConversion_list    = molidutils.uniquelist_to_conversionlist(unique_molID_occur_map)


#####################################
# INTERACTIVE MOLID CONVERSION MENU #
#####################################
# Goal: 
# users complete their molID conversion templates by ensuring that each member of
# molIDConversion_list has status complete = True
input_menu = ""
input_menu_complete = ""
# For use in the next section
concat_menu = ""

while(input_menu != "QUIT"):
    if (input_menu_complete == "1"):
        print("""Congratulations! You have successfully constructed your
                 conversion templates. You can proceed to the next section 
                 by selection option 7 or, continue to edit your conversion 
                 template through this menu
              """)
    print("""PDBClean MolID Conversion Build Menu
             Select one of the following options to proceed:
             1) Show full conversion
             2) Show only unassigned conversions
             3) Enter input file
             4) Search MolID to add chain ID conversion
             5) Go entry by entry to add chain ID conversion
             6) Remove a chain ID conversion
          """)
    if (input_menu_complete == "1"):
        print("    7) Continue to next step of curation")
    input_menu = input('Option Number: ')
    if (input_menu == "1"):
        molidutils.show_full_conversion(molIDConversion_list)
    elif (input_menu == "2"):
        molidutils.show_unassigned_conversion(molIDConversion_list)
    elif (input_menu == "3"):
        molIDConversion_list = molidutils.add_user_conversion(molIDConversion_list)
    elif (input_menu == "4"):
        molIDConversion_list = molidutils.edit_conversion_interface(molIDConversion_list, action='add')
    elif (input_menu == "5"):
        molIDConversion_list = molidutils.edit_conversion_manual(molIDConversion_list)
    elif (input_menu == "6"):
        molIDConversion_list = molidutils.edit_conversion_interface(molIDConversion_list, action='remove')
    elif (input_menu == "7"):
        if (input_menu_complete == "1"):
            input_menu = "QUIT"
            concat_menu = "START"
    input_menu_complete = molidutils.check_complete(molIDConversion_list)

########################################
# INTERACTIVE MOLID CONCATENATION MENU #
########################################
# Goal: 

if (concat_menu == "START"):
    # Prepare for concatenation step
    # We now have to take the information contained in the MolIDConversion objects
    # in molIDConversion_list to update the MolID objects in master_molID_class_list
    # We then need to mine these updated MolID objects to figure out which ones
    # contain concatenated chains. These will be presented to the user in another
    # interactive menu section where they can update the planned conversion on
    # a file by file basis
    
    master_molID_class_list = molidutils.update_masterlist(master_molID_class_list, molIDConversion_list)

    concat_menu = ""
    concat_menu_complete = ""

    while(concat_menu != "QUIT"):

        count_problems = molidutils.problem_counter(master_molID_class_list)
        if (count_problems == 0):
            concat_menu_complete = "1"

        if (concat_menu_complete == "1"):
            print("""Congratulations! You have successfully constructed your
                     conversion templates.You can proceed to the next section 
                     by selection option 7 or, continue to edit your conversion 
                     template through this menu
                  """)
        print("""PDBClean Concatenations Menu
                 Note: All proposed concatenations must be accepted before the curation can
                 be completed.
                 Select one of the following options to proceed:
                 1) Show all conversions
                 2) Show only unaccepted concatenations
                 3) Search and modify destination chainIDs of proposed concatenations
                 4) Search and modify order of proposed concatenations
                 5) Search and accept proposed concatenations
              """)
        if (concat_menu_complete == "1"):
            print("    6) Finalize Curation")

        concat_menu = input('Option Number: ')

        if (concat_menu == "1"):
            molidutils.show_full_conversion(master_molID_class_list, step='concatenation')
        elif (concat_menu == "2"):
            molidutils.show_unassigned_conversion(master_molID_class_list, step='concatenation')
        elif (concat_menu == "3"):
            master_molID_class_list = molidutils.edit_concatenation_interface(master_molID_class_list, action='try')[0]
        elif (concat_menu == "4"):
            master_molID_class_list, new_order = molidutils.edit_concatenation_interface(master_molID_class_list, action='update')
        elif (concat_menu == "5"):
            master_molID_class_list = molidutils.edit_concatenation_interface(master_molID_class_list, new_order=new_order, action='accept')[0]
        elif (concat_menu == "6"):
            print("Finalizing Curation ...")
            molidutils.masterlist_to_pdb(filelist, master_molID_class_list, target_dir=target_dir)
            concat_menu = "QUIT"
