import os, glob
import re
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

#
def process(projdir=None, step='clean', source='raw_bank', target='clean_bank', pdbformat='.cif', verbose=True):
    """
    process
    """
    if projdir is not None:
        source_dir = projdir+'/'+source
        target_dir = projdir+'/'+target
        input_list = glob.glob(source_dir+'/*'+pdbformat)
        i=0
        for input_cif in input_list:
            cif_name=os.path.basename(input_cif)
            if verbose:
                i+=1
                print('[{0}/{1}]: {2}'.format(i,len(input_list),cif_name))
            output_cif=target_dir+'/'+cif_name
            if(step=='clean'):
                if os.path.isfile(output_cif):
                    os.remove(output_cif)
                clean_cif(input_cif, output_cif) #FAPA, CLOSE )
            elif(step=='simplify'):
                # missing line: remove all assembly cif already created
                simplify_cif(input_cif, output_cif, pdbformat)
            elif(step=='fixhet'):
                fixhet_cif(input_cif, output_cif, pdbformat)
            elif(step=='finalize'):
                finalize(input_cif, output_cif, pdbformat)
#
def finalize(oldfile, newfile, pdbformat):
    """
    finalize
    """
    with open(oldfile) as myfile:
        newciffile = open(newfile, 'w')
        for line in myfile:
            line_split = line.strip()
            line_split = line.split()
            if (line_split[0] == "ATOM") or (line_split[0] == "HETATM"):
                newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[17] + " " + line_split[7] + " " + line_split[8] + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + line_split[15] + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
            else:
                newciffile.write(line)
#
def fixhet_cif(oldfile, newfile, pdbformat):
    """
    fixhet_cif
    """
    with open(oldfile) as myfile:
        hetch_count_map = {}
        newciffile = open(newfile, 'w')
        last_res = 0
        last_res = {}
        usage = {}
        for line in myfile:
            line_split = line.strip()
            line_split = line_split.split()
            # Can Update so that it looks at all atoms?
            if (line_split[0] == "HETATM"):
                key = str(line_split[15]) + "_" + str(line_split[16]) + "_" + str(line_split[17]) + "_" + str(line_split[18])
                if line_split[17] in hetch_count_map:
                    if (line_split[15] != last_res[line_split[17]]) or (key in usage):
                        if key not in usage:
                            usage[key] = 1
                        hetch_count_map[line_split[17]] += 1
                        newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
                        last_res[line_split[17]] = line_split[15]
                    else:
                        usage[key] = 1
                        newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
                else:
                    usage[key] = 1
                    hetch_count_map[line_split[17]] = 1
                    last_res[line_split[17]] = line_split[15]
                    newciffile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + " " + line_split[3] + " " + line_split[4] + " " + line_split[5] + " " + line_split[6] + " " + line_split[7] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[9] + " " + line_split[10] + " " + line_split[11] + " " + line_split[12] + " " + line_split[13] + " " + line_split[14] + " " + str(hetch_count_map[line_split[17]]) + " " + line_split[16] + " " + line_split[17] + " " + line_split[18] + " " + line_split[19] + "\n")
            else:
                newciffile.write(line)
#
def simplify_cif(oldfile, newfile, pdbformat):
    """
    simplify_cif
    """
    mmcif_dict = MMCIF2Dict(oldfile)
    #
    # Create map from asym_id to assembly_id
    # assembly_id information may be either str or list type, so I'm going to force
    # to list type
    asym_assembly_map = {}
    assembly_id = mmcif_dict['_pdbx_struct_assembly_gen.assembly_id']

    if not isinstance(assembly_id, list):
        assembly_id_list = []
        asym_id_list = []
        assembly_id_list.append(assembly_id)
        asym_id_list.append(mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list'])
    else:
        assembly_id_list = []
        assembly_id_list = assembly_id
        asym_id_list = mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list']
    # Convert asym_id entry into a list of asym_ids
    for i in range(len(assembly_id_list)):
        asym_id = asym_id_list[i]
        asym_id = asym_id.strip()
        asym_id = re.sub(' ', '', asym_id)
        asym_id = asym_id.split(',')
        for ident in asym_id:
            asym_assembly_map[ident] = assembly_id_list[i]

    for assembly in assembly_id_list:

        if (len(assembly_id_list)==1):
            newciffilename = str(re.sub(pdbformat, '', newfile))+"+00"
        else:
            newciffilename = str(re.sub(pdbformat, '', newfile))+"+0"+str(assembly)
        newciffile = open(newciffilename+pdbformat, 'w')
        newciffile.write("data_"+newciffilename+"\n")
        #Write entry.id
        newciffile.write("#\n")
        L = str(mmcif_dict['_entry.id']) # FAPA: Change to STR because the default was a list
        entryid = '_entry.id   ' + L
        newciffile.write(entryid + "\n")
        # Write Audit category
        newciffile.write("#\n")
        newciffile.write("loop_\n")
        newciffile.write("_citation_author.name\n")
        L = mmcif_dict['_citation_author.name']
        if isinstance(L, list):
            for i in L:
                newciffile.write("'" + re.sub("'", "", i) + "'" + "\n")
        else:
            newciffile.write("'" + re.sub("'", "", L) + "'" + "\n")
        # Write Citation category
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_citation.title" + "\n")
        newciffile.write("_citation.year" + "\n")
        newciffile.write("_citation.pdbx_database_id_DOI" + "\n")
        L1 = mmcif_dict['_citation.title']
        L2 = mmcif_dict['_citation.year']
        L3 = mmcif_dict['_citation.pdbx_database_id_DOI']
        if isinstance(L1, list):
            for i in range(len(L1)):
                newciffile.write("'" + re.sub("\n"," ",re.sub("'", "", L1[i])) + "' " + L2[i] + " " + L3[i] + "\n") #FAPA
        else:
            newciffile.write("'" + re.sub("\n"," ",re.sub("'", "", L1[i])) + "' " + L2[i] + " " + L3[i] + "\n") #FAPA
        # Write Resolution category
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_exptl.method" + "\n")
        newciffile.write("_exptl.resolution" + "\n")
        if '_exptl.method' in mmcif_dict:
            L1 = mmcif_dict['_exptl.method']
        elif '_refine_hist.pdbx_refine_id' in mmcif_dict:
            L1 = mmcif_dict['_refine_hist.pdbx_refine_id']
        else:
            L1 = mmcif_dict['_refine.pdbx_refine_id']
        if '_refine.ls_d_res_high' in mmcif_dict:
            L2 = mmcif_dict['_refine.ls_d_res_high']
        elif '_em_3d_reconstruction.resolution' in mmcif_dict:
            L2 = mmcif_dict['_em_3d_reconstruction.resolution']
        elif '_refine_hist.d_res_high' in mmcif_dict:
            L2 = mmcif_dict['_refine_hist.d_res_high']
        else:
            L2 = '????'
        if isinstance(L1, list) and isinstance(L2, list):
            for i in range(len(L1)):
                newciffile.write("'" + L1[i] + "' " + L2[i] + " " + "\n")
        elif isinstance(L1, list) and not isinstance(L2,list):
            newciffile.write("'" + L1[0] + "' " + L2 + " " + "\n")
        elif not isinstance(L1,list) and isinstance(L2,list):
            newciffile.write("'" + L1 + "' " + L2[0] + " " + "\n")
        else:
            newciffile.write("'" + L1 + "' " + L2 + " " + "\n")
        # Write Entity category
        # This was the part fixed by FAPA
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_entity.id" + "\n")
        newciffile.write("_entity.pdbx_description" + "\n")
        L1 = mmcif_dict['_entity.id']
        L2 = mmcif_dict['_entity.pdbx_description']
        for i in range(len(L1)):
            L2[i] = L2[i].upper()
            L2[i] = L2[i].replace(":", "")
            newciffile.write(L1[i] + " '" + L2[i].replace("'", "") + "'\n")

        # Now write the coordinate portion of the file
        newciffile.write("#" + "\n")
        newciffile.write("loop_" + "\n")
        newciffile.write("_atom_site.group_PDB" + "\n")
        newciffile.write("_atom_site.id" + "\n")
        newciffile.write("_atom_site.type_symbol" + "\n")
        newciffile.write("_atom_site.label_atom_id" + "\n")
        newciffile.write("_atom_site.label_alt_id" + "\n")
        newciffile.write("_atom_site.label_comp_id" + "\n")
        newciffile.write("_atom_site.label_asym_id" + "\n")
        newciffile.write("_atom_site.label_entity_id" + "\n")
        newciffile.write("_atom_site.label_seq_id" + "\n")
        newciffile.write("_atom_site.pdbx_PDB_ins_code" + "\n")
        newciffile.write("_atom_site.Cartn_x" + "\n")
        newciffile.write("_atom_site.Cartn_y" + "\n")
        newciffile.write("_atom_site.Cartn_z" + "\n")
        newciffile.write("_atom_site.occupancy" + "\n")
        newciffile.write("_atom_site.B_iso_or_equiv" + "\n")
        newciffile.write("_atom_site.auth_seq_id" + "\n")
        newciffile.write("_atom_site.auth_comp_id" + "\n")
        newciffile.write("_atom_site.auth_asym_id" + "\n")
        newciffile.write("_atom_site.auth_atom_id" + "\n")
        newciffile.write("_atom_site.pdbx_PDB_model_num" + "\n")
        L1 = mmcif_dict['_atom_site.group_PDB']
        L2 = mmcif_dict['_atom_site.id']
        L3 = mmcif_dict['_atom_site.type_symbol']
        L4 = mmcif_dict['_atom_site.label_atom_id']
        L5 = mmcif_dict['_atom_site.label_alt_id']
        L6 = mmcif_dict['_atom_site.label_comp_id']
        L7 = mmcif_dict['_atom_site.label_asym_id']
        L8 = mmcif_dict['_atom_site.label_entity_id']
        L9 = mmcif_dict['_atom_site.label_seq_id']
        L10 = mmcif_dict['_atom_site.pdbx_PDB_ins_code']
        L11 = mmcif_dict['_atom_site.Cartn_x']
        L12 = mmcif_dict['_atom_site.Cartn_y']
        L13 = mmcif_dict['_atom_site.Cartn_z']
        L14 = mmcif_dict['_atom_site.occupancy']
        L15 = mmcif_dict['_atom_site.B_iso_or_equiv']
        L16 = mmcif_dict['_atom_site.auth_seq_id']
        L17 = mmcif_dict['_atom_site.auth_comp_id']
        L18 = mmcif_dict['_atom_site.auth_asym_id']
        L19 = mmcif_dict['_atom_site.auth_atom_id']
        L20 = mmcif_dict['_atom_site.pdbx_PDB_model_num']
        # This is strictly to clean up a common problem found in ribosome structures
        # HETATMs are defined mid chain, ie in an asymmetric unit of non-HETATMs
        # Want to treat these atoms different from other HETATMs
        asymatom_list = []
        for i in range(len(L1)):
            if L1[i] == "ATOM":
                if L7[i] not in asymatom_list:
                    asymatom_list.append(L7[i])
        #
        hetatm_map = {}
        hetatm_count = 0

        ### BEGIN NEW CODE : THIS SECTION FIXES THE BUG
        ### This section is necessary to print the biological assemblies on separate files

        BioAssembly = mmcif_dict['_pdbx_struct_assembly_gen.asym_id_list']

        for i in range(len(L1)):
            if (L7[i] in list(re.sub(",", "", BioAssembly[int(assembly)-1]))): # Only print the chains pertaining to particular biological assembly unit
                if (L1[i]=="ATOM"):
                    newciffile.write(L1[i] + " " + L2[i] + " " + L3[i] + ' "' + L4[i] + '" ' + L5[i] + " " + L6[i] + " " + L7[i] + " " + L8[i] + " " + L9[i] + " " + L10[i] + " " + L11[i] + " " + L12[i] + " " + L13[i] + " " + L14[i] + " " + L15[i] + " " + L16[i] + " " + L17[i] + " " + L18[i] + ' "' + L19[i] + '" ' + L20[i] + "\n")
                elif (L1[i]=="HETATM"):
                    if L7[i] in asymatom_list:
                        newciffile.write("ATOM" + " " + L2[i] + " " + L3[i] + ' "' + L4[i] + '" ' + L5[i] + " " + L6[i] + " " + L7[i] + " " + L8[i] + " " + L9[i] + " " + L10[i] + " " + L11[i] + " " + L12[i] + " " + L13[i] + " " + L14[i] + " " + L15[i] + " " + L16[i] + " " + L17[i] + " " + L18[i] + ' "' + L19[i] + '" ' + L20[i] + "\n")
                    else:
                        newciffile.write(L1[i] + " " + L2[i] + " " + L3[i] + ' "' + L4[i] + '" ' + L5[i] + " " + L6[i] + " " + L7[i] + " " + L8[i] + " " + L9[i] + " " + L10[i] + " " + L11[i] + " " + L12[i] + " " + L13[i] + " " + L14[i] + " " + L15[i] + " " + L16[i] + " " + L17[i] + " " + "het" + L6[i] + L18[i] + ' "' + L19[i] + '" ' + L20[i] + "\n")
        newciffile.write("#" + "\n")


        ### END NEW CODE

#
def clean_cif(oldfile, newfile):
    """
    clean_cif
    """
    entry_list = ['_entry.id',
                  '_atom_site.group_PDB',
                  '_citation_author.name',
                  '_citation.title',
                  '_pdbx_struct_assembly_gen.assembly_id',
                  '_entity.pdbx_description',
                  '_exptl.method',
                  '_em_3d_reconstruction.resolution',
                  '_refine_hist.pdbx_refine_id',
                  '_refine.pdbx_refine_id']
    keylength_list = [ 9,
                       20,
                       21,
                       15,
                       37,
                       24,
                       13,
                       32,
                       27,
                       22]
    with open(oldfile) as old_file:
        alllines = []
        linecount = 0
        poundline = 0
        flag = 0
        for line in old_file:
            alllines.append(line)
            if linecount == 0:
                with open(newfile, 'a') as new_file:
                    new_file.write(alllines[0])
            for entry, keylength in zip(entry_list, keylength_list):
                flag = check_and_write_entry(entry, line, alllines, line[0:keylength], flag, range(poundline, linecount), newfile)
            if '#' in line[0]:
                poundline = linecount
            linecount += 1
        with open(newfile, 'a') as new_file:
            new_file.write('#\n')
#
def check_and_write_entry(entry, line, alllines, key, flag, linerange, newfile):
    """
    check_and_write_entry
    """
    if entry in key:
        flag = 1
    elif (flag==1) and '#' in line[0]:
        with open(newfile, 'a') as new_file:
            for i in linerange:
                new_file.write(alllines[i])
        flag=0
    return flag
