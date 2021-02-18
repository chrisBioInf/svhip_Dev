# Svhip - fast and flexible ncRNA-classifier training

## Introduction

Svhip is a software developed in Python for the automatic retraining of the RNAz 2.0 SVM-Classifier, a bioinformatics tool for the de novo discovery of non-coding RNA in whole genome screens.   
It processes Clustal-Alignments or raw FASTA files and returns a list of feature tuples, that can either be used to recalibrate or test the RNAz software or train a completely independent classifier. Working examples for the latter are included in the repository (https://github.com/chrisBioInf/svhip_Dev), provided is a Random Forests, a Logistic Regression classifier and others.  
The RNAz classification software can be downloaded here: https://www.tbi.univie.ac.at/software/RNAz/#download.  

The basic Svhip-workflow looks something like this:  

    Alignment/FASTA file 
         |                       |<---add more data <--|<----recalibration <------|           
         |                       |                     |                          |
         |                       |                     |                          |
         --------> Feature vector creation ----> Classifier Training ------> Testing 
                                                                                  |
                                                                                  |
                                                                                  |
                                                                          Use Classifer on
                                                                             genomic data


## Installation

Follow the steps outlined below to install the current version of Svhip. First, clone the git repository: 

```bash
git clone https://github.com/chrisBioInf/svhip_Dev.git
```
Make sure you have a current version of Python3 and pip installed:
```bash
python --version
python3 -m pip --version
```
If not already installed, you will need the `wheel` utility for installation:
```bash
python -m pip install --upgrade pip setuptools wheel
```

Navigate to the cloned git repository:
```bash
cd Svhip_Dev
```
And install the current Svhip distribution with:
```bash 
python3 -m pip install svhip_dev.whl
```

{TODO: Integrate with python-application 2.8.0 for simplicity of use. This weekend.}

## Basic Usage

To run Svhip, simply type: 
```bash
svhip
```

Since no other arguments are given, this will direct you to the manual page. Take note of the following basic command line arguments: 
```text
svhip -i [INPUT] -o [OUTPUT] [OPTION] [OPTION ARGUMENTS]

where viable OPTION parameters are:

	-testset: 	Takes the arguments given in [OPTION ARGUMENTS] as well as an 
			[INPUT] parameter and writes a .dat file of training feature vectors
			based on input files. Please note that testset_creation module
			has its own set of options, for a complete list refer to help 
			function of this module, i.e.:
			>> svhip -testset --help
			In any case, alignment windows are created and can be later used as
			testing instance with RNAz.
	-write_m:	Takes the arguments given in [OPTION ARGUMENTS] as well as an 
			[INPUT] parameter, optimizes hyper parameters and writes a .model 
			file based on training data
			gathered using the testset_creation module (a .dat file written by
			the -testset or -auto methods). Please note that the
			write_model module has its own set of options, for a complete 
			list refer to help function of this module, i.e.:
			>> svhip -write_m --help
	-auto:		Combines both above options. Takes an [INPUT] parameter and creates 
			a .model file. Note that at the moment this option will be
			executed using the default parameters for both modules, so 
			the input syntax amounts to:
			>> svhip -i [INPUT] -o [OUTPUT FILE] -auto
			Alignment windows created are still saved for potential later use.

```
These three form the core functionality of Svhip. `testset` will take an alignment file you provide with the `-i [alignment name]` parameter and calculate a set of corresponding feature vectors for training. These will by default be written to a file called `[alignment name].dat` but can be specified with `-o [output name]`, like:
```bash
svhip -i [alignment] -o myoutput.dat -data_gen
```
Svhip will automatically determine the input format (clustal or FASTA file).
In general, Svhip utilizes four main suffixes to designate different file types used during the process:
```text
	.dat		.dat files contain calculated features of alignments/alignment windows. 
			They are created any time a new input alignment is processed and may contain
			positive and negative or only positive feature vector instances. These are 
			already scaled, can be directly used for classifier training using the 
			-write_m command line argument and are utilized to hold training data 
			for quick later recalibration or to allow flexible addition of new data.
			For this reason they are also saved as an intermediate result when following
			the full pipeline from alignment to trained classifier using the -auto command.
			Manually editing them is highly deprecated. 

	.model		Contain all data necessary for the trained classifier to function and can be simply
			loaded using the edited RNAz software. Loading a new .model file with RNAz is
			equivalent to swapping he (dinucleotide, sequence based) classifier. End product of
			the pipeline.

	.aln_align_n_i	Partial alignments, where n refers to the number of sequences per alignment and
			i is an internal index counting the number of files created. These contain the alignment 
			windows created from the pre processed input alignment (default columns: 120, step size: 40).
			From these, feature vectors contained in the .dat files are calculated. They also serve 
			as test set instances, as they are already processed for prediction using RNAz 
			(see RNAz documentation).
	
	.random		'Negative' control alignments, created with SISSIz v. 0.11. These serve as dinucleotide and 
			gap-controlled null models for both negative instance creation and structural conservation
			validation. Remember to also include the randomized alignment windows when constructing a 
			test set with Svhip. For classifier training purposes, they are automatically included 
			upon generation. 

```

Note that these are primarily used to distinguish file types and can be renamed if the need arises.   
With that in mind, let's return to the usage example. If we want to use the newly generated feature vector file (aln.dat from now on) to train an external classifier, we are already done: We can export the plain text file and use it in any way we want (examples for a python-based parsing is provided in the repository along example alternate classifiers). If we want to use it to train a classifier for RNAz 2.0, we would type:

```bash
svhip -i aln.dat -o aln.model -write_m 
``` 

This will suffice to write all necessary classifier data in the `aln.model` file once training is finished. This can then immediatly be used with the RNAz version edited for dynamic model loading or by replacing the native dinucleotide_decision.model file before compilation. Please also have a look at the individual `--help` pages of the `testset` and `write_m` arguments that further explain individual options, like maximum pairwise sequence identity, depth of crossvalidation during classifier training and others.  
If you want to create a `.model` file for RNAz and want to take a shortcut, you can also combine both command explained above, by typing:

```bash
svhip -i [alignmentname] -o [modelfilename] -auto
```
This will effectively combine both commands and immediatly generate a functional `.model` file. All function parameters from both `write_m` and `data_gen` can be used with `auto`. 


## Contact

Christopher Klapproth  
christopher@bioinf.uni-leipzig.de  
University Leipzig   
Interdisciplinary Centre for Bioinformatics  

## License
[MIT](https://choosealicense.com/licenses/mit/)
