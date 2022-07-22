from __future__ import print_function
from __future__ import division
import os
import time

# AA Map from 3 letter amino acid id to 1 letter id
def ResnConvert(resn):
    AA = {}
    AA["UNK"] = "X"
    AA["ALA"] = "A"
    AA["ARG"] = "R"
    AA["ASN"] = "N"
    AA["ASP"] = "D"
    AA["CYS"] = "C"
    AA["GLN"] = "Q"
    AA["GLU"] = "E"
    AA["GLY"] = "G"
    AA["HIS"] = "H"
    AA["ILE"] = "I"
    AA["LEU"] = "L"
    AA["LYS"] = "K"
    AA["MET"] = "M"
    AA["PHE"] = "F"
    AA["PRO"] = "P"
    AA["SER"] = "S"
    AA["THR"] = "T"
    AA["TRP"] = "W"
    AA["TYR"] = "Y"
    AA["VAL"] = "V"
    AA["A"]   = "A"
    AA["C"]   = "C"
    AA["U"]   = "U"
    AA["G"]   = "G"
    if resn not in AA:
        ans = "X"
    else:
        ans = AA[resn]
    return ans
# END AA Map from 3 letter amino acid id to 1 letter id

def AlignSequences(sequence_vec):
    # Takes a list of sequence strings and performs a MUSCLE alignment, outputting a vector of aligned sequence strings
    with open("Seq.fa", 'w') as newfafile:
        for seq in sequence_vec:
            newfafile.write("> Seq" + "\n")
            newfafile.write(seq + "\n")

    os.popen('muscle -in Seq.fa -out Seq.afa')

    time.sleep(1)

    aligned_seq = []
    with open("Seq.afa") as seqfile:
        seq = ""
        for line in seqfile:
            if (line[0] == ">"):
                if (seq != ""):
                    aligned_seq.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        aligned_seq.append(seq)
    return (aligned_seq)
# END AlignSequences

def AlignSequences_v2(sequence_vec, file_name):
    # Takes a list of sequence strings and performs a MUSCLE alignment, outputting a vector of aligned sequence strings
    with open(file_name+".fa", 'w') as newfafile:
        i = 0
        for seq in sequence_vec:
            newfafile.write("> Seq " + str(i)+ "\n")
            newfafile.write(seq + "\n")
            i += 1
    command = "muscle -in "+file_name+".fa -out "+file_name+".afa"
    os.popen(command)
    time.sleep(1)
    aligned_seq_map = {}
    aligned_seq = []
    seq = ""
    with open(file_name + ".afa") as seqfile:
        for line in seqfile:
            if (line[0] == ">"):
                # Very first line
                if (seq == ""):
                    line = line.strip()
                    line = line.split()
                    key = line[2]
                else:
                    aligned_seq_map[key] = seq
                    seq = ""
                    line = line.strip()
                    line = line.split()
                    key = line[2]
            else:
                seq += line.strip()
        aligned_seq_map[key] = seq
    for i in range(len(aligned_seq_map)):
        aligned_seq.append(aligned_seq_map[str(i)])
    return (aligned_seq)
# END AlignSequences



def ScoreSequenceAlignment(seq1, seq2):
    # Scores based on exact identity. Should  maybe be updated to take longer
    # of sequences so that it can be used with unaligned seq strings too
    score = 0
    for i in range(len(seq1)):
        if (seq1[i] == seq2[i]):
            score += 1
    score = score/len(seq1)
    return score

