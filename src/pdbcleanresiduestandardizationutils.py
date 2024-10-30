from __future__ import print_function
from __future__ import division
from Bio.PDB.MMCIFParser import FastMMCIFParser
from PDBClean.alignmentutils import *
from PDBClean.listutils import *
#

####################
# INITIALIZE STEPS #
####################

def pdb_to_structurelists(filelist):
    """
    Iterates through a list CIF(s) and retrieves structure IDs, chain IDs, and maps chain IDs to their sequences,
    and maps chain IDs to their residue numbers.

    Parameters:
    -----------
    filelist : list of str
    	list of file paths for all '.cif' files in specified directory

    Returns:
    --------
    Structure_Sequences : dict
        Contains dictionary where chain ID is mapped to their sequence for each structure.
    ChID_ResiNum_Vector : list of dict
        Each dictionary maps the chain ID to their residue numbers for a structure
    structid_list : list of str
    	List of unique structure identifiers for each CIF. Format is 'input directory / CIF'
    chid_list : list of str
    	A list containing all the chain IDs from CIF(s)
    """

    # Structure_Sequences is the master list of maps from chain IDs to their sequences
    Structure_Sequences = {}
    ChID_ResiNum_Vector = []
    chid_list = []
    structid_list = []

    N = 0
    for my_file in filelist:
        N += 1
        print("Reading:" + ' ' + my_file + "  (" + str(N) + " of " + str(len(filelist)) + ")")
        struct = FastMMCIFParser(auth_residues=False,QUIET=1).get_structure(str(my_file), my_file)
        structid_list.append(struct.get_id())
        chid_seq_map = {}
        chid_resinum_map = {}
        # Only written for structures with only one model in them
        for chain in struct[0]:
            if (chain.get_id() not in chid_resinum_map): #FAPA: HERE WE NEED TO ADD IF TO CHECK IF pdbx_PDB_ins_code != '?'
                chid_resinum_map[chain.get_id()] = []
            key = str(struct.get_id()) + "_" + str(chain.get_id())
            resinum_list = []
            seq = ""
            for residue in chain:
                #print("printing residue ids: "+ str(residue.get_id()[2]))
                resinum_list.append(residue.get_id()[1])
                #chid_resinum_map[chain.get_id()].append(residue.get_id()[1])
                # For each residue we extract both the residue number and the associated "letter" (pdbx_PDB_ins_code)
                chid_resinum_map[chain.get_id()].append(str(residue.get_id()[1])+str(residue.get_id()[2])) #FAPA 17 oct 2024
                #resinum_list.append(residue.id()[1])
                #chid_resinum_map[chain.get_id()].append(residue.id()[1])
                seq += ResnConvert(residue.get_resname())
            Structure_Sequences[key] = seq
            # chid_list is a master list of all chainIDs used
            chid_list.append(chain.get_id())
        ChID_ResiNum_Vector.append(chid_resinum_map)
    chid_set = set(chid_list)
    chid_list = sorted(list(chid_set))
    #print(ChID_ResiNum_Vector)
    return Structure_Sequences, ChID_ResiNum_Vector, structid_list, chid_list

#########################################
# INTERACTIVE STANDARDIZATION FUNCTIONS #
#########################################

def perform_multiple_alignment(Structure_Sequences, ChID_ResiNum_Vector, structid_list, chid_list, check):
    """
    Interactive user interface for performing multiple alignments

    Parameters:
    -----------
    Structure_Sequences : dict
        Contains dictionary where chain ID is mapped to their sequence for each structure.
    ChID_ResiNum_Vector : list of dict
        Each dictionary maps the chain ID to their residue numbers for a structure
    structid_list : list of str
    	List of unique structure identifiers for each CIF. Format is 'input directory / CIF'
    chid_list : list of str
    	A list containing all the chain IDs from CIF(s)
    check : str
        Option chosen by user which opens the submenu

    Returns:
    --------
    Structure_Sequences_Aligned : dict
        A dictionary where each key is a combination of structure identifier and chain ID,
        and the value is the aligned sequence for that chain.
    Structure_ConversionTemplate : dict
        A dictionary mapping each structure identifier to a conversion template,
        which contains mappings of residue numbers from the original sequence to the aligned sequence.
    chid_list : list of str
        Updated list of chain IDs where some may have been removed based on the user's options
    check : str
        Updated string representing the state of the main menu, set to '1' to indicate a state change.
    """
    Structure_Sequences_Aligned = {}
    Structure_ConversionTemplate = {}
    Structure_Sequences_GAPS = {}
    input_submenu = ""
    while(input_submenu != "QUIT"):
        print("    Perform multiple alignments to identify residues",
              "    1) Show list of chains to be standardized",
              "    2) Remove chain IDs from list of chains to be standardized",
              "    3) Input file of chain IDs to remove from list of chains to be standardized",
              "    4) Perform multiple alignments",
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
                this_chainsseq_list = []
                this_chainsseq_list_ids = [] #FAPA
                this_chainsseq_aligned_list = []
                for I in range(len(structid_list)):
                    key = str(structid_list[I]) + "_" + chid
                    if key in Structure_Sequences:
                        this_chainsseq_list.append(Structure_Sequences[key])
                        this_chainsseq_list_ids.append(structid_list[I]) # FAPA
                #this_chainsseq_aligned_list_map = AlignSequences_v2(this_chainsseq_list, chid,this_chainsseq_list_ids )
                #this_chainsseq_aligned_list_map = AlignSequences_v3(this_chainsseq_list, chid, this_chainsseq_list_ids) #FAPA MAY2024
                this_chainsseq_aligned_list_map, this_chainseq_gap_percentages = AlignSequences_v4(this_chainsseq_list, chid,
                                                                    this_chainsseq_list_ids)  # FAPA JULY2024
                i = 0
                for I in range(len(structid_list)):
                    key = str(structid_list[I]) + "_" + chid
                    if key in Structure_Sequences:
                        #Structure_Sequences_Aligned[key] = this_chainsseq_aligned_list[i]
                        Structure_Sequences_Aligned[key] = this_chainsseq_aligned_list_map[str(structid_list[I])]
                        Structure_Sequences_GAPS[key] = this_chainseq_gap_percentages #FAPA
                        i += 1

            #THIS IS THE VERSION THAT WORKS, COMMENTED SO WE TRY SOMETHING NEW
            #for I in range(len(structid_list)):
            #    conversion_template = {}
            #    for chain in ChID_ResiNum_Vector[I]:
            #        resinum_aligned_list = []
            #        key = str(structid_list[I]) + "_" + str(chain)
            #        if key in Structure_Sequences_Aligned:
            #            seq = Structure_Sequences_Aligned[key]
            #            i = 0
            #            for resn in seq:
            #                i += 1
            #                if (resn != "-"):
            #                    resinum_aligned_list.append(i)
            #            i = 0
            #            for residue in ChID_ResiNum_Vector[I][chain]:
            #                key2 = chain + "_" + str(residue)
            #                conversion_template[key2] = resinum_aligned_list[i]
            #                i += 1
            #    Structure_ConversionTemplate[structid_list[I]] = conversion_template

            # MY TEST STARTS HERE, WITH VARIATIONS OF THE CODE ABOVE

            for I in range(len(structid_list)):
                conversion_template = {}
                for chain in ChID_ResiNum_Vector[I]:
                    #print(ChID_ResiNum_Vector[I])
                    #residue_numbers_users = [residue.get_id()[1] for residue in chain.get_residues()]
                    #print(residue_numbers_users)
                    resinum_aligned_list = []
                    key = str(structid_list[I]) + "_" + str(chain)
                    if key in Structure_Sequences_Aligned:
                        seq = Structure_Sequences_Aligned[key]
                        gaps = Structure_Sequences_GAPS[key]
                        i = 0
                        counter=1
                        new_res_num=[]
                        freq_tracker=1
                        gap_tracker=0
                        gap_letter=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
                        for freq in gaps:
                            #freq_tracker = 1
                            #gap_tracker = 0
                            #print(freq)
                            if freq < 30:
                                new_res_num.append(counter)
                                counter += 1
                                freq_tracker=freq
                                gap_tracker=0
                                #print("hello")
                            else:
                                #print(gap_tracker)
                                #new_res_num.append(str(counter)+"_"+str(gap_tracker))
                                new_res_num.append(str(counter-1) +" "+str(gap_letter[gap_tracker]))
                                gap_tracker+=1
                                freq_tracker=freq

                        #print(new_res_num)
                        #print(len(new_res_num))
                        #print(len(seq))


                        for resn in seq:
                            #print(new_res_num[i])
                            if (resn != "-"):
                                resinum_aligned_list.append(new_res_num[i])
                            i += 1

                        #print("resinum_aligned_list")
                        #print(len(resinum_aligned_list))
                        #print(resinum_aligned_list)

                        #i = 0

                        #print("THIS IS THE CHAIN NUMBER TEST")
                        #print(ChID_ResiNum_Vector[I][chain]) #FAPA

                        for residue in range(len(resinum_aligned_list)):
                            #i +=1
                            #print("this is the value of i:", i)
                            #print("this is resimnum_aligned_list[i]:", resinum_aligned_list[residue] )
                            #print("this is the value of residue:", residue)
                            #print("this is what is in the structure:",ChID_ResiNum_Vector[I][chain][residue])

                            key2 = chain + "_" + str(ChID_ResiNum_Vector[I][chain][residue])
                            #print("key 2 is:",key2)
                                #print(key2) #FAPA
                                #print(resinum_aligned_list[i]) #FAPA
                            conversion_template[key2] = resinum_aligned_list[residue]
                            i += 1
                        #print(conversion_template)



                        #for residue in ChID_ResiNum_Vector[I][chain]:
                        #    #i +=1
                        #    print("this is the value of i:", i)
                        #    print("this is resimnum_aligned_list[i]:", resinum_aligned_list[i] )
                        #    print("this is the value of residue:", residue)
                        #    if i == len(resinum_aligned_list)-1:
                        #        break
                        #    elif resinum_aligned_list[i] in ChID_ResiNum_Vector[I][chain]:
                        #        key2 = chain + "_" + str(residue)
                        #        #print(key2) #FAPA
                        #        #print(resinum_aligned_list[i]) #FAPA
                        #        conversion_template[key2] = resinum_aligned_list[i]
                        #    i += 1
                        #print(conversion_template)

                Structure_ConversionTemplate[structid_list[I]] = conversion_template




            check = "1"
            input_submenu = "QUIT"
    return Structure_Sequences_Aligned, Structure_ConversionTemplate, chid_list, check

def show_conversiontemplate(Structure_ConversionTemplate):
    """
    Prints the conversion template to screen

    Paramters:
    ----------
    Structure_ConversionTemplate : dict
        A dictionary mapping each structure identifier to a conversion template,
        which contains mappings of residue numbers from the original sequence to the aligned sequence.

    Returns:
    --------
    None
    """

    for structid in Structure_ConversionTemplate:
        print(structid)
        for key in Structure_ConversionTemplate[structid]:
            print(key + ":" + str(Structure_ConversionTemplate[structid][key]))

## FAPA
def write_and_show_conversiontemplate(Structure_ConversionTemplate, target_dir, write_csv=True):
    """
    Writes and displays a mapping of old residue IDs to new residue IDs for each structure.

    This function prints the mapping to the console and optionally writes it to a CSV file
    in the specified target directory.

    Parameters:
    -----------
    Structure_ConversionTemplate : dict
    	list of file paths for all '.cif' files in specified directory
    target_dir : str
        Directory where the new files will be saved
    write_csv : bool, optional
        Writes the mapping to a CSV file named 'OldResID_NewResID_Map.csv' if True. Default is 'True'.

    Returns:
    --------
    None
    """

    if write_csv:
        with open(f'{target_dir}/OldResID_NewResID_Map.csv', 'w') as fout:
            fout.write('OldResID:NewResId:File\n')


    for structid in Structure_ConversionTemplate:
        #print(structid)
        for key in Structure_ConversionTemplate[structid]:
            #print(key + ":" + str(Structure_ConversionTemplate[structid][key]))
            #structid_for_print=[x.split("/")[-1] for x in structid]
            structid_for_print = structid.split("/")[-1]
            if write_csv:
                with open(f'{target_dir}/OldResID_NewResID_Map.csv', 'a') as fout:
                    fout.write(f'{key}:{str(Structure_ConversionTemplate[structid][key])}:{structid_for_print}\n')

# FAPA

#################
# FINALIZE STEP #
#################

def conversiontemplate_to_pdb(filelist, Structure_ConversionTemplate, target_dir=None):
    """
    Saves the conversion template into re-written CIF(s) which are placed into the target directory

    Parameters:
    -----------
    filelist : str
    	list of file paths for all '.cif' files in specified directory
     Structure_ConversionTemplate : dict
        A dictionary mapping each structure identifier to a conversion template,
        which contains mappings of residue numbers from the original sequence to the aligned sequence.
    target_dir : str, optional
        Directory where the new files will be saved. If none, no files will be saved.

    Returns:
    --------
    None
    """
    for my_files in filelist:
        newciffilename=target_dir+'/'+my_files.split('/')[-1]
        with open(my_files) as myfile:
            with open(newciffilename, 'w') as newciffile:
                # Now figure out which file is which template
                conversion_template = Structure_ConversionTemplate[myfile.name]
                for line in myfile:
                    if (line[0:4] == "ATOM") or (line[0:6] == "HETATM"):
                        # Chains outside map should not exist but just in case
                        line_split = line.strip()
                        line_split = line.split()
                        key = line_split[17] + "_" + str(line_split[15])
                        if key in conversion_template:
                            newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(conversion_template[key]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n"
                            newciffile.write(newline)
                        else:
                            newciffile.write(line)
                    else:
                        newciffile.write(line)

# FAPA MAY 2024 TEST STARTS

def conversiontemplate_to_pdb_FAPA(filelist, Structure_ConversionTemplate, target_dir=None):
    """
    Saves the conversion template into re-written CIF(s) which are placed into the target directory.
    This function considers cases where a residue number also includes a letter.

    Parameters:
    -----------
    filelist : str
    	list of file paths for all '.cif' files in specified directory
     Structure_ConversionTemplate : dict
        A dictionary mapping each structure identifier to a conversion template,
        which contains mappings of residue numbers from the original sequence to the aligned sequence.
    target_dir : str, optional
        Directory where the new files will be saved. If none, no files will be saved.

    Returns:
    --------
    None
    """
    for my_files in filelist:
        newciffilename=target_dir+'/'+my_files.split('/')[-1]
        with open(my_files) as myfile:
            with open(newciffilename, 'w') as newciffile:
                # Now figure out which file is which template
                conversion_template = Structure_ConversionTemplate[myfile.name]
                #print(conversion_template)


                for line in myfile:
                    resnum=1
                    #old_line_resnum=0
                    if (line[0:4] == "ATOM") or (line[0:6] == "HETATM"):
                        # Chains outside map should not exist but just in case
                        #line_split = line.strip()
                        line_split = line.split()
                        #print(line_split[8])
                        #new_line_resnum=int(line_split[8])

                        #if new_line_resnum == start_line_resnum:

                        # We need to consider the value of pdbx_PDB_ins_code, in column 9
                        # This is considered in the key
                        # and original value will be overwritten with '?'
                        # in next version, we will add the letter value to column 9

                        if str(line_split[9]) == '?':
                            key = line_split[6] + "_" + str(line_split[8]) + " "  # FAPA: WE WANT CHAINID_RESID TO BE THE KEY
                        else:
                            key = line_split[6] + "_" + str(line_split[8]) + str(line_split[9])


                        #key = line_split[6] + "_" + str(line_split[8]) #FAPA: WE WANT CHAINID_RESID TO BE THE KEY
                        #key = line_split[6] + "_" + str(resnum)  # FAPA: WE WANT CHAINID_RESID TO BE THE KEY

                        #print(len(line_split))



                        if key in conversion_template:
                            #print(key, conversion_template[key])
                            if len(line_split) == 18:
                                if len(str(conversion_template[key]).split()) < 2:
                                    #print(conversion_template[key])
                                    newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + \
                                              line_split[7] + " " + str(conversion_template[key]) + " " + "?" + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + \
                                              " " + line_split[14] + " " + str(conversion_template[key]) + " " + line_split[16] + " " + line_split[17] + " " + "\n"
                                    newciffile.write(newline)
                                else:
                                    # the key already contains two columns (this is a test)
                                    #print("i am here: "+ str(conversion_template[key].split()[1]))
                                    newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + \
                                              line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + \
                                              line_split[6] + " " + \
                                              line_split[7] + " " + str(conversion_template[key].split()[0]) + " " + \
                                              str(conversion_template[key].split()[1])+ " " + \
                                              line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + \
                                              line_split[13] + \
                                              " " + line_split[14] + " " + str(conversion_template[key].split()[0]) + " " + \
                                              line_split[16] + " " + line_split[17] + " " + "\n"
                                    newciffile.write(newline)
                            else:
                                if len(str(conversion_template[key]).split()) < 2:
                                    #print(conversion_template[key])
                                    newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + \
                                              line_split[7] + " " + str(conversion_template[key]) + " " + "?" + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + \
                                              " " + line_split[14] + " " + str(conversion_template[key]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n"
                                    newciffile.write(newline)
                                else:
                                    newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + \
                                              line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + \
                                              line_split[6] + " " + \
                                              line_split[7] + " " + str(conversion_template[key].split()[0]) + " " + \
                                              str(conversion_template[key].split()[1]) + " " + \
                                              line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + \
                                              line_split[13] + \
                                              " " + line_split[14] + " " + str(conversion_template[key].split()[0]) + " " + \
                                              line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n"
                                    newciffile.write(newline)
                        else:
                            if len(line_split) == 18:
                                newciffile.write(line)
                            else:
                                if line_split[8] == ".":
                                    newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + \
                                              line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + \
                                              line_split[6] + " " + \
                                              line_split[7] + " " + line_split[15] + " " + \
                                              line_split[9] + " " + \
                                              line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + \
                                              line_split[13] + \
                                              " " + line_split[14] + " " + line_split[15] + " " + \
                                              line_split[16] + " " + line_split[17] + " " + line_split[18] + " "+ line_split[19] +"\n"
                                    newciffile.write(newline)
                                else:
                                    newciffile.write(line)
                    else:
                        newciffile.write(line)

# FAPA MAY 2024 TEST ENDS