# PDBClean-0.0.2

PDBClean allows you to create a curated ensemble of structures deposited in the Protein Data Bank.

# What does PDBClean actually do?

We have created Jupyter Notebooks to show how to create a curated ensemble of structures using PDBClean.

### [Step 0. Download structural ensemble form RCSB PDB.](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step0.DownloadStructuralEnsembleFromRCSBPDB.ipynb)

Download all structures that match the name and sequence of your molecule of interest.

> **Note:** This notebook does not display on the Github website, download and open in your browser.   

### [Step 1. Clean Structures and Create one CIF file per biological assembly.](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step1.CreateOneCIFFilePerBiologicalAssembly.ipynb)

A CIF file may contain multiple biological assemblies within one asymmetric unit. In this step we separate these biological assemblies, and create one CIF file for each one. We also reduce the number of data blocks included in the CIF file.

### [Step 2.1. Assign MOLID to the entities found in the CIF files, version 1](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step2.1.AssignMolIDToEntitiesFoindInCIFfiles1.ipynb)

The script goes over all the CIF files and collects all entities. The user can decide what Mol ID to assign them. In this example, we show the case in which we give a different ID to each entity found.
This step is also important because it lists all the entities that were found in your ensemble, so it allows you to identify if there is a structure that doesn't belong. We show an example of this in this notebook.

### [Step 2.2. Assign MOLID to the entities found in the CIF files, version 2](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step2.2.AssignMolIDToEntitiesFoindInCIFfiles2.ipynb)

Same as Step 2.1, but in our example, we give the same MOL ID to different entities. You may want to do this for example, if you want to give the same MOL ID to all ligands, or water molecules. Doing this will trigger a concatenation menu, which we show how to use.

### [Step 3. Chain ID standardization](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step3.ChainIDStandardization.ipynb)

Step 2 allows us to name each entity with whatever name we want. Step 3 makes sure that the chains that are the same (we do sequence alignment to determine similarity) in different cif files, have a consistent name.

### [Step 4. Residue ID Standardization](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/Step4.ResidueIDStandardization.ipynb)

Following step 3, now that we have consistent chain (entity) naming among all structures in the ensembe, we want to make sure that the numbering is also consistent (that the same residue position has the same number in all structures).

This is also the last step! You have a curated dataset!


> **Note:** There are more advanced curation steps and analysis that we will cover in future notebooks.

### [Check project mini tutorial](https://github.com/fatipardo/PDBClean-0.0.2/blob/master/Notebooks/CheckProject_CheckCreateDelete.ipynb)

This mini tutorial can be run after doing step 2. It shows you how our check_project tool works.

# Installation

We recommend installing PDBClean inside a virtual environment. We provide an `environment.yml` with the libraries you will need. Additionally, Anaconda is a recommended prerequisite before utilizing PDBClean. You will need to link muscle to your virtual environment, just follow the instructions we provide.
We have tested the installation on MacOS.

1. Download PDBClean from GitHub and install environment from YML file

>git clone git@github.com:fatipardo/PDBClean-0.0.2.git

>cd PDBClean-0.0.2

>conda env create -f environment.yml

2. Activate environment and install PDBClean

>conda activate PDBCleanV2

>python setup.py install

3. Download Muscle and Link Muscle to our new environment

PDBClean uses [Muscle](https://drive5.com/muscle5/)
(click the link and download the newest version for your computer). In order to use Muscle, we must first change the file's permissions, to allow it to execute.

>chmod a+x muscle5.1.{YOUR VERSION}
- Example: chmod a+x muscle5.1.macos_intel64
- If you are using MacOS, you will get a warning, click `Cancel`. Then go to your computers `System Preferences`. Open `Security & Privacy`. Then `General`. You will see a notification about muscle, click `Allow Anyway`. 

Now you should be able to run Muscle by typing in the terminal:
>./muscle5.1.{YOUR VERSION}

In order to link Muscle to your virtual environment, type in the terminal:

>ln -s {PATH/TO/MUSCLE}/muscle  {PATH/TO/ANACONDA}/anaconda3/envs/PDBCleanV2/bin/
- Example: ln -s ~/Downloads/muscle5.1.macos_intel64 /Users/username/opt/anaconda3/envs/PDBCleanV2/bin/muscle

4. Install Jupyter Notebook kernel

python -m ipykernel install --user --name PDBCleanV2 --display-name PDBCleanV2


5. Running notebook:

> cd Notebooks

> jupyter notebook

- Open any notebook you would like to run.
- If Jupyter does not recognize the kernel, select ‘PDBCleanV2’ from the drop down menu.


## PDBClean team

The code in this repository is based on the code found [here](https://test.pypi.org/project/PDBClean/#files).
The code was originally written by Frédéric Poitevin and Nicholas Corsepius.
Fátima Pardo Avila and Liv Weiner created this repository.
We all worked on this project while being part of the Levitt Lab at Stanford University.
