from __future__ import print_function
from __future__ import division
import sys
import argparse
import re
import numpy as np
import csv
import os
import copy
import time
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from PDBClean.alignmentutils import *
from PDBClean.listutils import *

####################
# INITIALIZE STEPS #
####################

def pdb_to_structurelists(filelist):
    """
    pdb_to_structurelists
    """
    Structure_Sequences = []
    Standard_Sequences = {}
    chid_list = []
    structid_list = []
    N = 0
    for my_file in filelist:
        N += 1
        print("Reading:" + ' ' + my_file + "  (" + str(N) + " of " + str(len(filelist)) + ")")
        struct = FastMMCIFParser(QUIET=1).get_structure(str(my_file), my_file)
        structid_list.append(struct.get_id())
        chid_seq_map = {}
        for chain in struct[0]:
            seq = ""
            for residue in chain:
                seq += ResnConvert(residue.get_resname())
            seq = re.sub('X', '', seq)
            if (len(seq) > 4):
                chid_seq_map[chain.get_id()] = seq
                chid_list.append(chain.get_id())
        Structure_Sequences.append(chid_seq_map)
    chid_set = set(chid_list)
    chid_list = sorted(list(chid_set))
    return Structure_Sequences, structid_list, chid_list


#########################################
# INTERACTIVE STANDARDIZATION FUNCTIONS #
#########################################

def select_standard_seq_from_reference(Structure_Sequences, Standard_Sequences, structid_list, input_menu_check_1):
    """
    select_standard_seq_from_reference
    """
    input_submenu = ""
    input_submenu_check_1 = ""
    while(input_submenu != "QUIT"):
        print("    Select Standard Sequences from input structure",
              "    1) Show list of input structures",
              "    2) Select input structure",
              sep="\n")
        if (input_submenu_check_1 == "1"):
            print("    3) Remove chains from Standard Sequences",
                  "    4) Inspect chains in Standard Sequences",
                  "    5) Return to main menu",
                  sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            show_list(structid_list)
        elif (input_submenu == "2"):
            Standard_Sequences, input_menu_check_1, input_submenu_check_1 = select_input_structure(Structure_Sequences, 
                                                                                                   structid_list, 
                                                                                                   input_menu_check_1,
                                                                                                   input_submenu_check_1)
        elif (input_submenu == "3" and input_submenu_check_1 == "1"):
            Standard_Sequences = remove_chains_from_standard(Standard_Sequence)
        elif (input_submenu == "4" and input_submenu_check_1 == "1"):
            inspect_chains_in_standard(Standard_Sequences)
        elif (input_submenu == "5" and input_submenu_check_1 == "1"):
            input_submenu = "QUIT"
    return Standard_Sequences, input_menu_check_1

def create_standard_seq_from_consensus(Structure_Sequences, Standard_Sequences, chid_list, input_menu_check_1):
    """
    create_standard_seq_from_consensus
    """
    input_submenu = ""
    input_submenu_check_1 = ""
    while(input_submenu != "QUIT"):
        print("    Create Standard Sequences from consensus of input structures.",
              "    Type QUIT to return to the main menu.",
              "    1) Show list of chain IDs for Standard Sequences",
              "    2) Enter chain IDs to remove from list",
              "    3) Input file with list of chain IDs to remove",
              "    4) Create Standard Sequences from consensus of input structures",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            show_list(chid_list)
        elif (input_submenu == "2"):
            chid_list = remove_user_defined_chain_from_list(chid_list)
        elif (input_submenu == "3"):
            chid_list = remove_file_defined_chain_from_list(chid_list)
        elif (input_submenu == "4"):
            for chid in chid_list:
                Standard_Sequences = assign_standard_from_consensus(Structure_Sequences, Standard_Sequences, chid)
            input_submenu = "QUIT"
            input_menu_check_1 = "1"
    return Standard_Sequences, input_menu_check_1

def review_standard_seq(Structure_Sequences, Standard_Sequences):
    """
    review_standard_seq
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Review Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of chain IDs in Standard Sequences",
              "    2) Enter chain ID of Standard Sequence to inspect/edit",
              "    3) Enter chain ID and inspect consensus of all matching chains",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            show_list(Standard_Sequences)
        elif (input_submenu == "2"):
            input_submenu = input('Chain ID: ')
            if input_submenu in Standard_Sequences:
                print(Standard_Sequences[input_submenu])
        elif (input_submenu == "3"):
            chid= input('Chain ID: ')
            this_chainsseq_list, this_chainsseq_score = get_this_chainsseq_list(Structure_Sequences, chid, verbose=True)

def align_to_standard_seq(Structure_Sequences, Standard_Sequences, structid_list):
    """
    align_to_standard_seq
    """
    ignore_chid = []
    input_submenu= ""
    input_submenu_4_check_1 = ""
    while(input_submenu != "QUIT"):
        print("    Perform pairwise alignments against Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of structure chain IDs to ignore when pairwise aligning to the Standard Sequences",
              "    2) Enter chain IDs to add to ignore list",
              "    3) Input file with list of chain IDs to add to ignore list",
              "    4) Perform pairwise alignments against Standard Sequences and create conversion template",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            for chid in sorted(ignore_chid):
                print(chid)
        elif (input_submenu == "2"):
            print("    Enter chain ID to add to ignore list. When complete, enter DONE.")
            while(input_submenu != "DONE"):
                input_submenu = input('Chain ID: ')
                if (input_submenu != "DONE"):
                    ignore_chid.append(input_submenu)
        elif (input_submenu == "3"):
            input_submenu = input('File: ')
            if (os.path.isfile(input_submenu) == True):
                my_file = open(input_submenu)
                for line in my_file:
                    ignore_chid.append(line.strip())
            else:
                print("File does not exist.")
        elif (input_submenu == "4"):
            ChainReassignmentMapping_List = []
            ChainReassignmentScores_List = []
            for chid_seq_map in Structure_Sequences:
                Struct_ChainReassignmentMap = {}
                Struct_ChainReassignmentScores = {}
                std_chid_list = []
                chidlist_list = []
                scorelist_list = []
                for std_chid in Standard_Sequences:
                    std_chid_list.append(std_chid)
                    structchid_score_map = {}
                    # First, check to see if chid is already correctly assigned
                    if std_chid in chid_seq_map:
                        aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[std_chid]])
                        score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                        # If score is perfect, don't bother aligning all chains
                        if (score == 1):
                            structchid_score_map[std_chid] = score
                            for struct_chid in chid_seq_map:
                                if (struct_chid != std_chid):
                                    structchid_score_map[struct_chid] = 0
                        # If score is not perfect, align all chains
                        else:
                            for struct_chid in chid_seq_map:
                                if struct_chid not in ignore_chid:
                                    aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                    score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                    structchid_score_map[struct_chid] = score
                    # Chid of standard does not match any structure chid, so perform all alignments
                    else:
                        for struct_chid in chid_seq_map:
                            if struct_chid not in ignore_chid:
                                aligned_seq = AlignSequences([Standard_Sequences[std_chid], chid_seq_map[struct_chid]])
                                score = ScoreSequenceAlignment(aligned_seq[0], aligned_seq[1])
                                structchid_score_map[struct_chid] = score
                    scorelist_list.append(list(reversed(sorted(structchid_score_map.values()))))
                    chidlist_list.append(list(reversed(sorted(structchid_score_map, key=structchid_score_map.get))))
                # Relevant information to create chain reassingment map has now been gathered from each structure via alignment
                # std_chid_list, chidlist_list and scorelist_list have the same indexes
                # At a given location i, std_chid_list[i] is the standard_sequence chid,
                # chidlist_list[i] is the list of chain ids in the given structure aligned to that standard sequence
                # and, scorelist_list[i] is the list of corresponding scores of that alignment
                # Note that there is no reason to keep all scores of all chids instead of just the top
                # scoring one, however, it's done because it seems likely that information may
                # be valuable and worth incorporating later on
                # Construct the map, making sure to avoid conflicts by comparing scores
                # These maps are old to new
                # CHANGES MADE HERE TO MAKE SURE CONFLICTS DO NOT ARRISE BETWEEN UNASSIGNED CHAINIDs and
                for i in range(len(std_chid_list)):
                    if chidlist_list[i][0] not in Struct_ChainReassignmentMap:
                        Struct_ChainReassignmentMap[chidlist_list[i][0]] = std_chid_list[i]
                        Struct_ChainReassignmentScores[chidlist_list[i][0]] = scorelist_list[i][0]
                    else:
                        # I < i
                        I = std_chid_list.index(Struct_ChainReassignmentMap[chidlist_list[i][0]])
                        # Compare score of original and new entry. If new entry has better score, replace
                        if(scorelist_list[i][0] > scorelist_list[I][0]):
                            Struct_ChainReassignmentMap[chidlist_list[i][0]] = std_chid_list[i]
                            Struct_ChainReassignmentScores[chidlist_list[i][0]] = scorelist_list[i][0]
                            if scorelist_list[I][0] in Struct_ChainReassignmentScores:
                                del Struct_ChainReassignmentScores[chidlist_list[I][0]]
                                del Struct_ChainReassignmentMap[chidlist_list[I][0]]
                # Now for final check, make sure no conflicts by checking unused chid
                # Loop over all chid in the structure
                unused_chid = []
                for chid in chid_seq_map:
                    if chid not in Struct_ChainReassignmentMap:
                        unused_chid.append(chid)
                        # chid is unused . . . will not change. If it is in destination list, then there will be conflict
                for chid in unused_chid:
                    used_dest_chid = list(Struct_ChainReassignmentMap.values()) #FAPA
                    if chid in used_dest_chid:
                        # chid is unused and a destination so a conflict exists and needs to be given destination
                        II = Structure_Sequences.index(chid_seq_map)
                        print("Conflict in structure " + structid_list[II])
                        for old_chid in Struct_ChainReassignmentMap:
                            if (Struct_ChainReassignmentMap[old_chid] == chid):
                                while (input_submenu != "PASS"):
                                    print("Original chain " + str(old_chid) + " is being reassigned to new chain " + str(chid))
                                    print("This creates a conflict because original chain " + str(chid) + " has not been reassigned.")
                                    print("Where should original chain " + str(chid) + " be reassigned?")
                                    input_submenu = input('New chain ID: ')
                                    if (input_submenu not in unused_chid) and (input_submenu not in Struct_ChainReassignmentMap.values()):
                                        # input is a valid destination
                                        input_submenu = "PASS"
                                    else:
                                        print("Not a valid destination chain ID")
                        if (len(input_submenu) <= 2):
                            Struct_ChainReassignmentMap[chid] = input_submenu
                            Struct_ChainReassignmentScores[chid] = 0
                            used_dest_chid.append(input_submenu)
                        else:
                            Struct_ChainReassignmentMap[chid] = input_submenu[0:2]
                            Struct_ChainReassignmentScores[chid] = 0
                            used_dest_chid.append(input_submenu[0:2])
                ChainReassignmentMapping_List.append(Struct_ChainReassignmentMap)
                ChainReassignmentScores_List.append(Struct_ChainReassignmentScores)
                input_menu_check_2 = "1"
                input_submenu = "QUIT"
    return ChainReassignmentMapping_List, ChainReassignmentScores_List, input_menu_check_2

def assign_standard_from_consensus(Structure_Sequences, Standard_Sequences, chid):
    """
    """
    this_chainsseq_list, this_chainsseq_score = get_this_chainsseq_list(Structure_Sequences, chid)
    # Now you have an ordered by length list of sequences pertaining to particular chID_list
    # You want to grab the longest one with the higest score and save it as this chains standard sequences
    this_score_vec = []
    for seq in this_chainsseq_list:
        this_score_vec.append(this_chainsseq_score[seq])
    this_score_vec = sorted(this_score_vec)
    topscore = this_score_vec[len(this_score_vec)-1]
    for seq in this_chainsseq_list:
        if (topscore == this_chainsseq_score[seq]):
            Standard_Sequences[chid] = seq
            break
    return Standard_Sequences

def get_this_chainsseq_list(Structure_Sequences, chid, verbose=False):
    """
    """
    this_chainsseq_list = []
    # Here, we're going through each structure in Structure_Sequences, which holds chid to seq maps
    for chid_seq_map in Structure_Sequences:
        if chid in chid_seq_map:
            this_chainsseq_list.append(chid_seq_map[chid])
    # Now we have a vector containing all the sequences pertaining to chid we want
    # Below is CRAP!
    # this_chainsseq_list is now sorted by length.
    # Start at longest, see if next longest is containing, if so add
    this_chainsseq_score = {}
    # Initialize scores
    this_chainsseq_set = set(this_chainsseq_list)
    for seq in this_chainsseq_set:
        this_chainsseq_score[seq] = 0
    # For each unique seq, count number of times in list
    for seq in this_chainsseq_list:
        this_chainsseq_score[seq] += 1
    # Reduce list to uniq entries
    this_chainsseq_list = sorted(list(this_chainsseq_set))
    this_chainsseq_list.sort(key=len)
    this_chainsseq_list.reverse()
    # Determine containedness of a sequence within a longer sequence
    for i in range(len(this_chainsseq_list)):
        for j in range(len(this_chainsseq_list)-i-1):
            if(this_chainsseq_list[i+1+j] in this_chainsseq_list[i]):
                this_chainsseq_score[this_chainsseq_list[i]] += this_chainsseq_score[this_chainsseq_list[i+1+j]]
    if verbose:
        for seq in this_chainsseq_list:
            print("Number of structures: " + str(this_chainsseq_score[seq]))
            print(seq)
    return this_chainsseq_list, this_chainsseq_score

def inspect_chains_in_standard(Standard_Sequences):
    """
    inspect_chains_in_standard
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Inspect chains in Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of chain IDs in Standard Sequences",
              "    2) Enter chain IDs inspect",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            for chid in Standard_Sequences:
                print(chid)
        elif (input_submenu == "2"):
            input_submenu = input('Chain ID: ')
            if input_submenu in Standard_Sequences:
                print(Standard_Sequences[input_submenu])
            else:
                print(input_submenu + " is not a chain ID in Standard Sequences")

def remove_chains_from_standard(Standard_Sequences):
    """
    remove_chains_from_standard
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Remove chains from Standard Sequences. Type QUIT to return to the main menu.",
              "    1) Show list of chain IDs in Standard Sequences",
              "    2) Enter chain IDs to remove from list",
              "    3) Input file with list of chain IDs to remove",
              sep="\n")
        input_submenu = input('Option Number: ')
        if (input_submenu == "1"):
            show_list(Standard_Sequences)
        elif (input_submenu == "2"):
            Standard_Sequences = remove_user_defined_chain_from_list(Standard_Sequences)
        elif (input_submenu== "3"):
            Standard_Sequences = remove_file_defined_chain_from_list(Standard_Sequences)
    return Standard_Sequences

def select_input_structure(Structure_Sequences, structid_list, input_menu_check_1, input_submenu_1_check_1):
    """
    select_input_structure
    """
    input_submenu = ""
    while(input_submenu != "QUIT"):
        input_submenu = input('File name: ')
        if input_submenu in structid_list:
            Standard_Sequences = Structure_Sequences[structid_list.index(input_submenu)]
            input_submenu = "QUIT"
            input_submenu_1_check_1 = str(1)
            input_menu_check_1 = str(1)
        else:
            print("""File name entered was not in list of files.""")
    return Standard_Sequences, input_menu_check_1, input_submenu_1_check_1


#################
# FINALIZE STEP #
#################

def reassignedmaps_to_pdb(filelist, ChainReassignmentMapping_List, structid_list, target_dir=None):
    """
    reassignmedmaps_to_pdb
    """
    for my_files in filelist:
        newciffilename=target_dir+'/'+my_files.split('/')[-1]
        with open(my_files) as myfile:
            with open(newciffilename, 'w') as newciffile:
                # Now figure out which file is which template
                    I = structid_list.index(myfile.name)
                    for line in myfile:
                        if (line[0:4] == "ATOM"):
                            # Chains outside map should not exist but just in case
                            line_split = line.strip()
                            line_split = line.split()
                            if line_split[17] in ChainReassignmentMapping_List[I]:
                                newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + line_split[15] + " " + line_split[16] + " " + ChainReassignmentMapping_List[I][line_split[17]] + " " + line_split[18] + " " + line_split[19] + "\n"
                                newciffile.write(newline)
                            else:
                                newciffile.write(line)
                        else:
                            newciffile.write(line)

def reassignedmaps_to_log(ChainReassignmentMapping_List, ChainReassignmentScores_List, structid_list, target_dir=None, verbose=True):
    """
    reassignedmaps_to_log
    """
    logfilename=target_dir+'/ChainStandardizationRecord.txt'
    with open(logfilename, 'w') as recordfile:
        for I in range(len(structid_list)):
            for chid in ChainReassignmentMapping_List[I]:
                recordfile.write(structid_list[I] + ":" + chid + ":" + ChainReassignmentMapping_List[I][chid] + ":" + str(ChainReassignmentScores_List[I][chid]) + "\n")
    if verbose:
        for i in range(len(ChainReassignmentMapping_List)):
            ChainReassignmentMap = ChainReassignmentMapping_List[i]
            ChainReassignmentScore = ChainReassignmentScores_List[i]
            print(structid_list[i])
            for oldchid in ChainReassignmentMap:
                print(oldchid + " " + ChainReassignmentMap[oldchid] + " " + str(ChainReassignmentScore[oldchid]))
#
