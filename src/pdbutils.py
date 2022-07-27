#
import os
import shutil
import re
import numpy as np
from urllib.request import urlopen
from contextlib import closing, suppress
#
def download_pdb_from_metadata(metadata, projdir=None):
    """
    download_pdb_from_metadata
    """
    if projdir is None:
        print("Please provide a project directory ...")
    else:
        download_dir=projdir+'/raw_bank'
        if not os.path.exists(download_dir):
            os.mkdir(download_dir)
        idset = get_idset_from_metadata(metadata)
        for pdbid in idset:
            download_pdb_from_id(pdbid, download_dir=download_dir)

def download_pdb_from_id(pdbid, pdbformat='.cif', download_dir=None):
    """
    download-pdb_from_id
    """
    download_url='https://files.rcsb.org/download/'
    if download_dir is None:
        print("Please provide a directory where to store downloaded files...")
    else:
        target = download_dir+'/'+pdbid+pdbformat
        source = download_url+pdbid.upper()+pdbformat
        download_from_url(source, target)

def get_idset_from_metadata(metadata):
    """
    get_idset_from_metadata
    """
    idlist = []
    for elt in metadata:
        idlist.append(elt[1:5])
    return sorted(set(idlist))

#
def retrieve_sequence_from_PDB(keyword, mode='sequence', update=True, seqfile=None):
    """
    retrieve_sequence_from_PDB: outputs a list of sequences fom seqfile (that needs to be updated if not input)

    Parameters:
    ===========
    - keyword: anything really
    - mode: 'sequence': whether keyword should match a sequence
            'metadata': or metadata
    - update: (True/False) whether the list needs to be updated
    - seqfile: if None, update needs to be True. The list retrieved online.
    """
    if update:
        with suppress(FileNotFoundError):
            os.remove(seqfile) # remove existing seqfile if any
        seqfile = retrieve_seqfile(seqfile=seqfile)
    metadata = []
    sequence = []
    with open(seqfile) as f:
        nextline=False
        prevline='#'
        for line in f:
            if nextline:
                sequence.append(line)
                nextline=False
            else:
                hit = re.findall(keyword, line, flags=re.I)
                if hit:
                    if(mode=='sequence'):
                        metadata.append(prevline)
                        sequence.append(line)
                    elif(mode=='metadata'):
                        metadata.append(line)
                        nextline = True
            prevline=line
    return np.atleast_1d(sequence), np.atleast_1d(metadata)
#
def retrieve_seqfile(seqfile=None):
    """
    retrieve_seqfile
    """
    sequrl='ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt'
    if seqfile is None:
        seqfile='seqfile.txt'
    download_from_url(sequrl, seqfile)
    return seqfile
#
def download_from_url(source, target):
    """
    download_from_url
    """
    with closing(urlopen(source)) as r:
        with open(target, 'wb') as f:
            shutil.copyfileobj(r,f)
    print('wrote {0} from {1}'.format(target, source))
