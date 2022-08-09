#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import sys
import re
import numpy as np
import csv
import os
import copy
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


####################
# INITIALIZE STEPS #
####################

def pdb_to_masterlist(filelist):
    """
    pdb_to_masterlist
    """
    master_molID_class_list = []
    N=0
    for my_files in filelist:
        N += 1
        print("Reading:"+' '+my_files+"  ("+str(N)+" of "+str(len(filelist))+")")
        with open(my_files) as myfile:
            this_molID_class = make_MolID_cif(myfile)
            master_molID_class_list.append(this_molID_class)
    return master_molID_class_list

def uniquelist_to_conversionlist(unique_map):
    """
    uniquelist_to_conversionlist
    """
    molIDConversion_list = []
    for key in unique_map:
        molIDConversion_list.append(initial_MolIDConversion(key, unique_map[key]))
    return molIDConversion_list

def update_masterlist(master_molID_class_list, molIDConversion_list):
    """
    update_masterlist
    """
    molIDConvert_molID_chID_map = {}
    for molIDConversion in molIDConversion_list:
        molIDConvert_molID_chID_map[molIDConversion.molID] = molIDConversion.chID_list
    # Update master_molID_class_list
    for molID_class in master_molID_class_list:
        for molID in molID_class.molID_chID:
            molID_class.add_chID_newchID_map(molID, molIDConvert_molID_chID_map[molID])
    # Determine which MolID objects now contain conflicts (concatenations)
    for molID_class in master_molID_class_list:
        molID_class.check_for_concatenations()
    return master_molID_class_list

# The class containing all the information about each file neccessary to build
# a conversion template
class MolID(object):
        # file_name = ""
        # chID_newchID_map = {}
        # molID_chID = {}
        # concat_order_map = {}
        # complete_order = {}

    def __init__(self, file_name, chID_newchID_map, molID_chID,
                 concat_order, complete_order):
        self.file_name = file_name
        self.chID_newchID_map = chID_newchID_map
        self.molID_chID = molID_chID
        self.concat_order = concat_order
        self.complete_order = complete_order

    def add_chID_newchID_map(self, molID, newchID_list):
        N = 0
        for oldchID in self.molID_chID[molID]:
            self.chID_newchID_map[oldchID] = newchID_list[N]
            N += 1

    def check_for_concatenations(self):
        self.concat_order = {}
        # Figure out which newchIDs appear more than once which will be
        # concatenated
        usage = {}
        duplicates = {}
        for chID in self.chID_newchID_map:
            if self.chID_newchID_map[chID] not in usage:
                usage[self.chID_newchID_map[chID]] = 1
            else:
                duplicates[self.chID_newchID_map[chID]] = 1
        # Create a map of oldchID to newchID for those being concatenated
        usage = {}
        for chID in self.chID_newchID_map:
            if self.chID_newchID_map[chID] in duplicates:
                if self.chID_newchID_map[chID] in usage:
                    usage[self.chID_newchID_map[chID]] += 1
                else:
                    usage[self.chID_newchID_map[chID]] = 1
                self.concat_order[chID] = usage[self.chID_newchID_map[chID]]
        # Update complete_order
        for chID in self.chID_newchID_map:
            if chID in self.concat_order:
                self.complete_order[chID] = False
            else:
                self.complete_order[chID] = True
    # END check_for_concatenations

    # This will force the complete_order[chID] to the input "complete" which is either True or False
    def force_complete_order(self, chID, complete):
        self.complete_order[chID] = complete

    def update_concat_order(self, chID, neworder):
        otheroldchID = ""
        if chID in self.concat_order:
            newchID = self.chID_newchID_map[chID]
            for oldchID in self.chID_newchID_map:
                if (self.chID_newchID_map[oldchID] == newchID) and (oldchID != chID) and (self.concat_order[oldchID] == neworder):
                    otheroldchID = oldchID
            if (otheroldchID != ""):
                self.concat_order[otheroldchID] = self.concat_order[chID]
                self.concat_order[chID] = neworder

# END MolID class


# make_molID is the function that will grab all relevant information from each
# file input
def make_MolID(myfile):

    file_name = myfile.name
    molID_list = []
    chID_list = []
    molID_chID = {}
    chID_newchID_map = {}
    concat_order = {}
    complete_order = {}
    complete = False
    for line in myfile:
        if "COMPND" in line:
            if "MOLECULE:" in line or "MOLID:" in line:
                line_list = line.rsplit(':')
                molID = line_list[1].strip()
                # Remove leading and trailing spaces
                molID = re.sub(' \{0,\}$|^ \{0,\}', '', molID)
                # Remove any special characters
                molID = re.sub(':|;|’|“|”', '', molID)
                # Add single leading and trailing space and bookmark with ;
                molID = re.sub('^', '; ', molID)
                molID = re.sub('$', ' ;', molID)
                molID_list.append(molID.upper())
            if "CHAIN:" in line:
                line_list = line.rsplit(':')
                chID = line_list[1].strip()
                # Remove all spaces and ;
                chID = re.sub(' |;', '', chID)
                chID = chID.split(',')
                chID_list.append(chID)
    for i in range(len(chID_list)):
        if molID_list[i] not in molID_chID:
            molID_chID[molID_list[i]] = chID_list[i]
        else:
            for chid in chID_list[i]:
                    molID_chID[molID_list[i]].append(chid)
    my_molID_class = MolID(file_name, chID_newchID_map,
                           molID_chID, concat_order, complete_order)
    my_molID_class.check_for_concatenations()
    return my_molID_class
# END make_MolID


# make_molID is the function that will grab all relevant information from each
# file input
def make_MolID_cif(myfile):

    file_name = myfile.name
    chID_newchID_map = {}
    molID_chID = {}
    concat_order = {}
    complete_order = {}
    complete = False

    mmcif_dict = MMCIF2Dict(myfile)

    # CIF files contain entity_id's which are used to link molID and chID
    # Need the entity_id and auth_asym_id correspondence
    # Put into entity_chIDlist_map
    entity_list = mmcif_dict['_atom_site.label_entity_id']
    chID_list = mmcif_dict['_atom_site.auth_asym_id']
    entity_chIDlist_map = {}

    for i in range(len(entity_list)):
        if entity_list[i] not in entity_chIDlist_map:
            entity_chIDlist_map[entity_list[i]] = [chID_list[i]]
        else:
            if chID_list[i] not in entity_chIDlist_map[entity_list[i]]:
                entity_chIDlist_map[entity_list[i]].append(chID_list[i])

    entity_list = mmcif_dict['_entity.id']
    molID_list = mmcif_dict['_entity.pdbx_description']
    molID_list = [i.upper() for i in molID_list]

    for i in range(len(entity_list)):
        # Need this step because some MolIDs may be present whose chains have been removed
        if entity_list[i] in entity_chIDlist_map:
            for chID in entity_chIDlist_map[entity_list[i]]:
                if molID_list[i] not in molID_chID:
                    molID_chID[molID_list[i]] = [chID]
                else:
                    molID_chID[molID_list[i]].append(chID)

    my_molID_class = MolID(file_name, chID_newchID_map,
                           molID_chID, concat_order, complete_order)
    my_molID_class.check_for_concatenations()
    return my_molID_class
# END make_MolID_cif


# MolIDConversion Class. The idea here is to create an object that
# will be modified by the user in the first step to map a MolID to a given
# chain ID with the correct number of occurances. This is the master list of
# molID's
class MolIDConversion(object):
    molID = ""
    chID_list = []
    occur = 0
    complete = bool

    def __init__(self, molID, chID_list, occur, complete):
        self.molID = molID
        self.chID_list = chID_list
        self.occur = occur
        self.complete = complete

    def check_for_completeness(self):
        if (len(self.chID_list) >= self.occur):
            self.complete = True
        else:
            self.complete = False

    def add_chID(self, chID):
        # if chID not in chID_list:
        #    self.chID_list.append(chID)
        self.chID_list.append(chID)

    def add_chID_list(self, chID_list):
        for chID in chID_list:
            # if chID not in self.chID_list:
            #    self.chID_list.append(chID)
            self.chID_list.append(chID)

    def remove_chID_list(self, chID_list):
        for chID in chID_list:
            if chID in self.chID_list:
                self.chID_list.remove(chID)


def initial_MolIDConversion(molID, occur):
    chID_list = []
    complete = False
    my_molID_conversion = MolIDConversion(molID, chID_list, occur, complete)
    return my_molID_conversion


def add_chID_list(MolIDConversion, chID_list):
    for chID in chID_list:
        if chID not in MolIDConversion.chID_list:
            MolIDConversion.chID_list.append(chID)
    return MolIDConversion
#

# Compile information from each file to create a master map of molID to max
# number of occurances
def CreateMasterUniqueMolIDMap(molID_class_list):
    unique_molID_map = {}
    for my_molID_class in molID_class_list:
        for molID in my_molID_class.molID_chID:
            if (unique_molID_map.get(molID) is None):
                unique_molID_map[molID] = \
                    len(my_molID_class.molID_chID[molID])
            elif (unique_molID_map.get(molID) <
                    len(my_molID_class.molID_chID[molID])):
                unique_molID_map[molID] = len(my_molID_class.molID_chID[molID])
    return unique_molID_map

#
#
# Read input file function
def read_input_file(input_cnv_file):
    user_molID_chID_map = {}
    user_molID_list_v = []
    user_chID_list_v = []
    # Read user input file and create user_molID_chID_map
    if (input_cnv_file != ""):
        if (os.path.isfile(input_cnv_file) is True):
            # Reading of input file
            with open(input_cnv_file) as mycnv:
                for line in mycnv:
                    line = line.strip()
                    my_line = line.split(':')
                    molID = my_line[0]
                    chID = my_line[1].split(',')
                    user_molID_list_v.append(molID)
                    user_chID_list_v.append(chID)
            # Mining of input file for information

            # for i in range(len(user_molID_list_v)):
            #     if (user_molID_chID_map.get(user_molID_list_v[i]) is None):
            #         user_molID_chID_map[user_molID_list_v[i]] = user_chID_list_v[i]
            #     elif (user_molID_chID_map.get(user_molID_list_v[i]) is not None):
            #         for this_chID in user_chID_list_v[i]:
            #             if this_chID not in user_molID_chID_map[user_molID_list_v[i]]:
            #                 user_molID_chID_map[user_molID_list_v[i]].append(this_chID)

            for molID_i, chID_i in zip(user_molID_list_v, user_chID_list_v):
                if molID_i not in user_molID_chID_map:
                    user_molID_chID_map[molID_i] = chID_i
                else:
                    for this_chID in chID_i:
                        if this_chID not in user_molID_chID_map[molID_i]:
                            user_molID_chID_map[molID_i].append(this_chID)

        return user_molID_chID_map
# End Read user input file function


####################################
# INTERACTIVE CONVERSION FUNCTIONS #
####################################

def show_full_conversion(current_list, step='conversion'):
    """
    show_full_conversion
    """
    if(step=='conversion'):
        for molIDConversion in current_list:
            molIDCon_chID_list_forPrint = re.sub('\[|\]| |\'', '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)
    elif(step=='concatenation'):
        for molID_class in current_list:
            for molID in molID_class.molID_chID:
                for chID in molID_class.molID_chID[molID]:
                    if chID in molID_class.concat_order:
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":"
                              + str(molID_class.concat_order[chID]))
                    else:
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":" +
                              str(0))

def show_unassigned_conversion(current_list, step='conversion'):
    """
    show_unassigned_conversion
    """
    if(step=='conversion'):
        for molIDConversion in current_list:
            molIDConversion.check_for_completeness()
            if molIDConversion.complete is False:
                molIDCon_chID_list_forPrint = re.sub('\[|\]| |\'', '', str(molIDConversion.chID_list))
                print(str(molIDConversion.occur)+":"+str(molIDConversion.molID)+":"+molIDCon_chID_list_forPrint)
    elif(step=='concatenation'):
        for molID_class in current_list:
            for molID in molID_class.molID_chID:
                for chID in molID_class.molID_chID[molID]:
                    if molID_class.complete_order[chID] is False:
                        print(molID_class.file_name + ":" + molID + ":"
                              + chID + ":" +
                              molID_class.chID_newchID_map[chID] + ":"
                              + str(molID_class.concat_order[chID]))

def add_user_conversion(molIDConversion_list):
    """
    add_user_conversion:
    Add contents of user's input file to your MolIDConversion class and check for completeness
    """
    input_cnv_file = input('Conversion File: ')
    user_molID_chID_map =  read_input_file(input_cnv_file)
    for molIDConversion in molIDConversion_list:
        # This is currently strict inclusion but perhaps should be except out
        # of convenience to the user
        for key in user_molID_chID_map:
            if (str(key)==str(molIDConversion.molID)):
                molIDConversion.add_chID_list(user_molID_chID_map[key])
        molIDConversion.check_for_completeness()
    # !! molIDConversion_list has been updated
    return molIDConversion_list

def edit_conversion_interface(moldIDConversion_list, action='add'):
    """
    edit_conversion_interface
    """
    search_term = input('MolID search term: ')
    molIDConversion_list, search_molIDConversion_list = molidutils.search_conversion(molIDConversion_list, search_term)
    input_submenu = 0
    while (input_submenu != "DONE"):
        print("    1) Further narrow down search results")
        if(action=='add'):
            print("    2) Add chain ID to conversion templates")
        elif(action=='remove'):
            print("    2) Remove chain ID from conversion templates")
        input_submenu = input('Option Number: ')
        if (input_submenu == "QUIT"):
            input_submenu = "DONE"
        if (input_submenu == "1"):
            search_term = input('MolID search term: ')
            molIDConversion_list, search_molIDConversion_list = molidutils.search_again_conversion(molIDConversion_list, search_molIDConversion_list, search_term)
        elif (input_submenu == "2"):
            if(action=='add'):
                print("""Enter new chain IDs, comma separated, no spaces""")
            elif(action=='remove'):
                print("""Enter chain ID to be removed, comma separated, no spaces""")
            chID_list = input('Chain IDs: ')
            if (chID_list == "") or (chID_list == "QUIT"):
                pass
            else:
                molIDConversion_list, search_molIDConversion_list = molidutils.edit_chain_conversion(molIDConversion_list, search_molIDConversion_list, chID_list, action=action)
            input_submenu = "DONE"
    return molIDConversion_list

def search_conversion(molIDConversion_list, search_term):
    """
    search_conversion
    """
    search_molIDConversion_list = []
    for molIDConversion in molIDConversion_list:
        if search_term in molIDConversion.molID:
            search_molIDConversion_list.append(molIDConversion)
            molIDCon_chID_list_forPrint = re.sub('\[|\]| |\'', '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)
    for molIDConversion in search_molIDConversion_list:
        molIDConversion_list.remove(molIDConversion)
    return molIDConversion_list, search_molIDConversion_list

def search_again_conversion(molIDConversion_list, search_molIDConversion_list, search_term):
    """
    search_again_conversion
    """
    new_search_molIDConversion_list = []
    for molIDConversion in search_molIDConversion_list:
        # If search term is found on search_molIDConversion_list
        # then add to new temporary list, if not, add back to
        # molIDConversion_list
        if search_term in molIDConversion.molID:
            new_search_molIDConversion_list.append(molIDConversion)
            molIDCon_chID_list_forPrint = re.sub('\[|\]| |\'', '', str(molIDConversion.chID_list))
            print(molIDConversion.molID+":"+molIDCon_chID_list_forPrint)
        else:
            molIDConversion_list.append(molIDConversion)
    search_molIDConversion_list = new_search_molIDConversion_list
    return molIDConversion_list, search_molIDConversion_list

def edit_chain_conversion(molIDConversion_list,search_molIDConversion_list, chID_list, action='add'):
    """
    edit_chain_conversion
    """
    chID_list = chID_list.split(',')
    # Performing chain addition to searched terms
    for molIDConversion in search_molIDConversion_list:
        if(action=='add'):
            molIDConversion.add_chID_list(chID_list)
        elif(action=='remove'):
            molIDConversion.remove_chID_list(chID_list)
    # Adding modified molIDConversions back to molIDConversion_list
    for molIDConversion in search_molIDConversion_list:
        molIDConversion.check_for_completeness()
        molIDConversion_list.append(molIDConversion)
    return molIDConversion_list, search_molIDConversion_list

def edit_conversion_manual(molIDConversion_list):
    """
    edit_conversion_manual
    """
    print("Enter chain IDs for each of the following MolID.")
    print("Comma separated, no spaces")
    for molIDConversion in molIDConversion_list:
        chID_list = input(molIDConversion.molID+":")
        if (chID_list == "") or (chID_list == "QUIT"):
            pass
        else:
            chID_list = chID_list.split(',')
        molIDConversion.add_chID_list(chID_list)
        molIDConversion.check_for_completeness()
    return molIDConversion_list

def check_complete(molIDConversion_list):
    """
    check_complete
    """
    num_unassigned = 0
    for molIDConversion in molIDConversion_list:
        if molIDConversion.complete is False:
            num_unassigned += 1
    if (num_unassigned > 0):
        input_menu_complete = "0"
    else:
        input_menu_complete = "1"
    return input_menu_complete

#######################################
# INTERACTIVE CONCATENATION FUNCTIONS #
#######################################

def problem_counter(master_molID_class_list):
    """
    problem_counter
    """
    count_problems = 0
    for molID_class in master_molID_class_list:
        for chID in molID_class.complete_order:
            if molID_class.complete_order[chID] is False:
                count_problems += 1
    return count_problems

def edit_concatenation_interface(master_molID_class_list, new_order=None, action='accept'):
    """
    edit_concatenation_interface
    """
    concat_submenu = 0
    while(concat_submenu != "QUIT"):
        print("I am here, but should not be 0_0")
        search_term, concat_submenu = get_search_term(concat_submenu)
        if concat_submenu != "QUIT":
            found_molID_class_chID_map, molID_class_been_copied = search_chains(master_molID_class_list, search_term)
        while (concat_submenu != "QUIT"):
            if(action=='try'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Update new chain ID
                      """)
            elif(action=='update'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Update concatenation order
                      """)
            elif(action=='accept'):
                print("""Select one of the following options to proceed:
                         1) Perform new search
                         2) Accept planned concatenation
                      """)
            concat_submenu = input('Option Number: ')
            if (concat_submenu == "1"):
                break
            elif (concat_submenu == "2"):
                if(action=='try'):
                    new_order = input('New Chain ID: ')
                elif(action=='update'):
                    new_order = input('New concatenation order: ')
                found_molID_class_chID_map = edit_chain_order(found_molID_class_chID_map, new_order, action=action)
                if(action=='try'):
                    print_conflicts(found_molID_class_chID_map)
                    print("Would you like to accept or deny these changes?")
                    while (concat_submenu != "QUIT"):
                        concat_submenu = input('Enter ACCEPT or DENY: ')
                        if (concat_submenu == "ACCEPT"):
                            print("yeah, yeah, we are here 1")
                            master_molID_class_list = accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                        concat_submenu = "QUIT"
                        concat_menu = ""
                else:
                    print("yeah, yeah, we are here else")
                    master_molID_class_list = accept_newchain(master_molID_class_list, found_molID_class_chID_map)
                    concat_submenu = "QUIT"
    print("we are now about to leave")
    print(master_molID_class_list)
    return master_molID_class_list, new_order

def get_search_term(value):
    """
    get_search_term
    """
    search_ok = ""
    while (search_ok != "OK"):
        print("Search format - File:MolID:OldChain:NewChain:ConcatOrder")
        search_term = input('Search ')
        if search_term == "QUIT":
            search_ok = "OK"
            value = "QUIT"
        search_term = search_term.split(":")
        if (len(search_term) > 4):
            search_ok = "OK"
    return search_term, value

def search_chains(master_molID_class_list, search_term):
    """
    search_chains
    """
    search_term_file = search_term[0]
    search_term_molID = search_term[1]
    search_term_oldchID = str(search_term[2])
    search_term_newchID = str(search_term[3])
    search_term_concatorder = search_term[4]

    found_molID_class_chID_map = {}
    molID_class_been_copied = {}
    # Perform search
    # Hits to the search have to be copied using deepcopy which will
    # create a distinct object that is an exact copy of another object
    # Modifications will be made to the copy to make sure they do not
    # introduce new errors
    for molID_class in master_molID_class_list:
        if (search_term_file in molID_class.file_name) or (search_term_file == ""):
            for molID in molID_class.molID_chID:
                if (search_term_molID in molID) or (search_term_molID == ""):
                    for chID in molID_class.molID_chID[molID]:
                        if (search_term_oldchID == chID) or (search_term_oldchID == ""):
                            newchID = molID_class.chID_newchID_map[chID]
                            if (search_term_newchID == newchID) or (search_term_newchID == ""):
                                if chID in molID_class.concat_order:
                                    this_concatord = molID_class.concat_order[chID]
                                else:
                                    this_concatord = 0
                                if (search_term_concatorder == str(this_concatord)) or (search_term_concatorder == ""):
                                    print(molID_class.file_name + ":" + molID + ":"
                                          + chID + ":" +
                                          molID_class.chID_newchID_map[chID] + ":"
                                          + str(this_concatord))
                                    if molID_class in molID_class_been_copied:
                                        found_molID_class_chID_map[molID_class_been_copied[molID_class]].append(chID)
                                    else:
                                        test_molID_class = copy.deepcopy(molID_class)
                                        molID_class_been_copied[molID_class] = test_molID_class
                                        found_molID_class_chID_map[test_molID_class] = [chID]
    return found_molID_class_chID_map, molID_class_been_copied

def edit_chain_order(found_map, newchID, action='try'):
    """
    edit_chain_order
    """
    for molID_class in found_map:
        for chID in found_map[molID_class]:
            if(action=='try'):
                molID_class.chID_newchID_map[chID] = newchID
                #HAVE TO FORCE EXISTING CONCATENATIONS BACK TO OK
                molID_class.check_for_concatenations()
            elif(action=='update'):
                molID_class.update_concat_order(chID, newchID)
            elif(action=='accept'):
                molID_class.force_complete_order(chID, True)
    return found_map

def print_conflicts(found_map):
    """
    print_conflicts
    """
    print("Updating this new chain ID will lead to the following conflicting assingments")
    was_printed = {}
    for molID_class in found_map:
        for molID in molID_class.molID_chID:
            for chID in molID_class.molID_chID[molID]:
                if molID_class.complete_order[chID] is False:
                    if chID not in was_printed:
                        was_printed[chID] = 1
                        print(molID_class.file_name + ":" + molID + ":" + chID + ":" + molID_class.chID_newchID_map[chID] + ":" + str(molID_class.concat_order[chID]))

def accept_newchain(masterlist, found_map):
    """
    accept_newchain
    """
    usage = {}
    for updated_molID_class in found_map:
        for molID_class in masterlist:
            if (molID_class.file_name == updated_molID_class.file_name):
                if updated_molID_class not in usage:
                    usage[updated_molID_class] = 1
                    masterlist.remove(molID_class)
                    masterlist.append(updated_molID_class)
                    
    return masterlist

#################
# FINALIZE STEP #
#################

def masterlist_to_pdb(filelist, masterlist, target_dir=None):
    """
    masterlist_to_pdb
    """
    for my_files in filelist:
        newciffilename=target_dir+'/'+my_files.split('/')[-1]
        with open(my_files) as myfile:
            with open(newciffilename, 'w') as newciffile:
                for molID_class in masterlist:
                    if (molID_class.file_name == myfile.name):
                        for line in myfile:
                            if (line[0:4] == "ATOM") or (line[0:6]=="HETATM"):
                                # Chains outside map should not exist but just in case
                                line_split = line.strip()
                                line_split = line.split()
                                if line_split[17] in molID_class.chID_newchID_map:
                                    # Residues have to be renumbered due to concatenations
                                    if line_split[17] in molID_class.concat_order:
                                        residue_offset = (molID_class.concat_order[line_split[17]] - 1) * 1000
                                        new_resinum = int(line_split[15]) + int(residue_offset)
                                        heresthenew_resinum = int(new_resinum)
                                        newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + molID_class.chID_newchID_map[line_split[17]] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(heresthenew_resinum) + " " + line_split[16] + " " + molID_class.chID_newchID_map[line_split[17]] + " " + line_split[18] + " " + line_split[19] + "\n"
                                        newciffile.write(newline)
                                    else:
                                        newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + molID_class.chID_newchID_map[line_split[17]] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + line_split[15] + " " + line_split[16] + " " + molID_class.chID_newchID_map[line_split[17]] + " " + line_split[18] + " " + line_split[19] + "\n"
                                        newciffile.write(newline)

                                else:
                                    newciffile.write(line)
                            # elif (line[0:6] == "COMPND"):
                            #     if "CHAIN:" in line:
                            #         newline = line[0:17] + molID_class.chID_newchID_map[line[17]] + line[18:]
                            #         newciffile.write(newline)
                            #     else:
                            #         newciffile.write(line)
                            else:
                                newciffile.write(line)

