# PDBClean-0.0.2

PDBClean allows you to create a curated ensemble of structures deposited in the Protein Data Bank. 

# What does PDBClean actually do? 

We have created Jupyter Notebooks to show how to create a curated ensemble of structures using PDBClean. 

## [Step 0. Download structural ensemble form RCSB PDB.](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step0.DownloadStructuralEnsembleFromRCSBPDB.ipynb)

Download all structures that match the name and sequence of your molecule of interested. 

## [Step 1. Clean Structures and Create one CIF file per biological assembly.](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step1.CreateOneCIFFilePerBiologicalAssembly.ipynb)

A CIF file may contain multiple biological assemblies within one asymmetric unit. In this step we separe these biologica assemblies, and create one CIF file for each one. We also reduce the number of data blocks included in the CIF file. 

## [Step 2.1. Assign MOLID to the entities found in the CIF files, version 1](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step2.1.AssignMolIDToEntitiesFoindInCIFfiles1.ipynb) 

The script foes over ll the CIF files and collects all entities. The user can decide what Mol ID to assign them. In this example we show the case in which we give a different ID to each entity found. 
This step is also important because it lists all the entities that were found in your ensemble, so it allows you to identify if there is a structure that doesn't belong. We show an example of this in this notebook. 

## [Step 2.2. Assign MOLID to the entities found in the CIF files, version 2](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step2.2.AssignMolIDToEntitiesFoindInCIFfiles2.ipynb) 

Same as Step 2.1, but in our example, we give the same MOL ID to different entities. You may want to do this for example, if you want to give the same MOL ID to all ligands, or water molecules. Doing this will trigger a concatenation menu, which we show how to use. 

## [Step 3. Chain ID standardization](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step3.ChainIDStandardization.ipynb) 

Step 2 allows us to name each entity with whatever name we want. Step 3 makes sure that the chains that are the same (we do sequence alignment to determine similarity) in different cif files, have a consistent name. 

## [Step 4. Residue ID Standardization](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step4.ResidueIDStandardization.ipynb)

Following step 3, now that we have consistent chain (entity) naming among all structures in the ensembe, we want to make sure that the numbering is also consistent (that the same residue position has the same number in all structures).

This is also the last step! You have a curated dataset! 


> **Note:** There are more advanced curation steps and analysis that we will cover in future notebooks. 

# Installation

We recommend installing PDBClean inside a virtual environment. We provide an `environment.yml` with the libraries you will need. 
PDBClean also uses [muscle](https://drive5.com/muscle5/), you will need to link muscle to your virtual environment, just follow the instructions we provide. 
We have tested the installation on MacOS. 

(**Note:** we are working on this section, come back soon for more details) 


### PDBClean team

The code in this repository is based on the code found [here](https://test.pypi.org/project/PDBClean/#files).
Which was originally written by Frédéric Poitevin and Nicholas Corsepius. 
Fátima Pardo Avila and Liv Weiner created this repository. 
We all worked on this project while being part of the Levitt Lab at Stanford University. 


