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
        struct = FastMMCIFParser(QUIET=1).get_structure(str(my_file), my_file)
        structid_list.append(struct.get_id())
        chid_resinum_map = {}
        # Only written for structures with only one model in them
        for chain in struct[0]:
            if (chain.get_id() not in chid_resinum_map):
                #Every chain ID you put inside its own list in the dictionary
                chid_resinum_map[chain.get_id()] = []
            key = str(struct.get_id()) + "_" + str(chain.get_id())
            resinum_list = []
            seq = ""
            for residue in chain:
                resinum_list.append(residue.get_id()[1])
                chid_resinum_map[chain.get_id()].append(residue.get_id()[1])
                seq += ResnConvert(residue.get_resname())
            Structure_Sequences[key] = seq
            # chid_list is a master list of all chainIDs used
            chid_list.append(chain.get_id())
        ChID_ResiNum_Vector.append(chid_resinum_map)
    chid_set = set(chid_list)
    chid_list = sorted(list(chid_set))
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
                for I in range(len(structid_list)):
                    key = str(structid_list[I]) + "_" + chid
                    if key in Structure_Sequences:
                        this_chainsseq_list.append(Structure_Sequences[key])
                        this_chainsseq_list_ids.append(structid_list[I]) # FAPA
                this_chainsseq_aligned_list_map = AlignSequences_v2(this_chainsseq_list, chid,this_chainsseq_list_ids)
                i = 0
                for I in range(len(structid_list)):
                    key = str(structid_list[I]) + "_" + chid
                    if key in Structure_Sequences:
                        Structure_Sequences_Aligned[key] = this_chainsseq_aligned_list_map[str(structid_list[I])]
                        i += 1
            for I in range(len(structid_list)):
                conversion_template = {}
                for chain in ChID_ResiNum_Vector[I]:
                    resinum_aligned_list = []
                    key = str(structid_list[I]) + "_" + str(chain)
                    if key in Structure_Sequences_Aligned:
                        seq = Structure_Sequences_Aligned[key]
                        i = 0
                        for resn in seq:
                            i += 1
                            if (resn != "-"):
                                resinum_aligned_list.append(i)
                        i = 0
                        for residue in ChID_ResiNum_Vector[I][chain]:
                            key2 = chain + "_" + str(residue)
                            conversion_template[key2] = resinum_aligned_list[i]
                            i += 1
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
        print(structid)
        for key in Structure_ConversionTemplate[structid]:
            print(key + ":" + str(Structure_ConversionTemplate[structid][key]))
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
                        line_split = line.split()
                        key = line_split[17] + "_" + str(line_split[15])
                        if key in conversion_template:
                            newline = line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(conversion_template[key]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n"
                            newciffile.write(newline)
                        else:
                            newciffile.write(line)
                    else:
                        newciffile.write(line)
