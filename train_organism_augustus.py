#!/usr/bin/env python

"""This script is for training augustus for a new species. 

parameters wanted: 
1. path to augustus
2. directory to save files
3. Species
4. Genbankfile
5. Number for percentage train and test set (80 = 80% of data as training and 20% as testset)

It consists of a few steps: 
1. Make a direction of the new species
2. Copy and edit all files from the 'standard' predefined files (generic folder) to your own organism
3. Creation of training and test set of the genbank file
4. Training augustus by optimising and editing the parameters (optimising and e-training)
5. Check the prediction accuracy
"""

import sys
import subprocess
import os.path
import os
from Bio import SeqIO
import random

def newDirectory(path_to_aug, species):
	"""Creation of new directory with the desired species

	path_to_aug - string with the correct path to the augustus directory
	species - the desired species
	"""
	speciesfolder = "{}/config/species/{}".format(path_to_aug, species)
	print "Create species directory"
	if os.path.isdir(speciesfolder):
		print "Directory already exists"
		pass
	else:
		print "Creation of species directory"
		cmd = "mkdir {}".format(speciesfolder)
		e = subprocess.check_call(cmd, shell=True)
	return speciesfolder

def copyGenericFiles(path_to_aug, species, path_to_species):
	"""Copy generic files to species directory
	
	path_to_aug - string with the correct path to the augustus directory
	species - the desired species
	path_to_species - path to the desired species
	"""
	print "Copy generic files to species folder ..."
	cmd_copy_files = "cp -a {}/config/species/generic/* {}".format(path_to_aug, path_to_species)
	if not os.listdir(path_to_species):
		print "Files are copied"
		e = subprocess.check_call(cmd_copy_files, shell=True)
	else:
		print "files already found in directory:", os.listdir(path_to_species)


def changeFilenames(speciesfolder, species):
	"""Change file names in species folder

	speciesfolder - path to the species
	species - name species (string) 
	"""
	for filename in os.listdir(speciesfolder):
		if filename.startswith("generic"):
			newname = filename.replace("generic", species)
			os.rename(os.path.join(speciesfolder, filename), os.path.join(speciesfolder, newname))

def changeParameterFile(speciesfolder, species):
	"""Change parameter file based on the species

	speciesfolder - path to the species
	species - name species (string)
	"""
	with open("{}/{}_parameters.cfg".format(speciesfolder, species), "r+") as inFile:
		for line in inFile:
			if "generic" in line:
				inFile.write(line.replace("generic", species))
			else:

				inFile.write(line)

def appendList(indexlist, l):
	"""Appending to index list
	
	indexlist - indexlist
	l - position
	"""
	appendedList = []
	for i in indexlist:
		appendedList.append(l[i])
	return appendedList

def splitGenbankInTrainTest(gb, p):
	"""Split genbank file into training and testset

	gb - genbankfile
	p - percentage in trainingset
	"""
	print "Creation of training and test set"
	entries = list(SeqIO.parse(gb, 'genbank'))
	num_entries = len(entries)
	percentage = (len(entries) * p) / 100
	train = set(random.sample(list(set(range(len(entries)))), percentage))
	test = set(range(len(entries))) - train
	trainingset = appendList(list(train), entries)
	testset = appendList(list(test), entries)
	
	print "DONE"
	return trainingset, testset

def writeGB(path, train, test):
	"""Writing train and testset into directory

	train - trainingset
	test - testset
	"""
	print "Creation of output directory"
	if os.path.isdir("{}/output".format(path)):
                print "Directory already exists"
                pass
        else:
                cmd = "mkdir {}/output".format(path)
                e = subprocess.check_call(cmd, shell=True)
		print "DONE"
	
	print "Creation of training and testset files"

	if os.path.isfile("{}/output/trainingSet.gb".format(path)):
		print "Trainingfile already exists"
	else:
		print "Creation trainingfile"
		with open(os.path.join("{}/output".format(path), "trainingSet.gb"), "w") as handle:
	                SeqIO.write(train, handle, 'genbank')
		print "DONE"

	if os.path.isfile("{}/output/testSet.gb".format(path)):
		print "Testfile already exists"
	else:
		print "Creation testfile"
		with open(os.path.join("{}/output".format(path), "testSet.gb"), "w") as handle:
			SeqIO.write(test, handle, 'genbank')
		print "DONE"

def optimiseParameters(path_to_aug, species, path_to_training):
	"""Optimising augustus parameters

	path_to_aug - path to augustus
	species - name of species
	train - trainingsfile
	"""

	training = {}
	cmd = "{}/scripts/optimize_augustus.pl --species={} --metapars={}/config/species/{}/{}_metapars.cfg --aug_exec_dir={}/bin --AUGUSTUS_CONFIG_PATH={}/config {}/output/trainingSet.gb" \
	.format(path_to_aug, species, path_to_aug, species, species, \
	path_to_aug, path_to_aug, path_to_training)
	e = subprocess.check_call(cmd, shell=True)

def etraining(path_to_aug, species, path_to_save):
	"""Perform etraining

	path_to_aug - path to augustus
	species - speciesname
	path_to_save
	"""
	cmd = "{}/bin/etraining --species={} --AUGUSTUS_CONFIG_PATH={}/config {}/output/trainingSet.gb"\
	.format(path_to_aug, species, path_to_aug, path_to_save)
	e = subprocess.check_call(cmd, shell=True)

def checkAccuracy():
	"""Check the accuracy of your prediction by using your testfile

	
	"""
	cmd = "{}/bin/augustus --species={} {}/output/testSet.gb"\
	.format(path_to_aug, species, testfile)
	e = subprocess.check_call(cmd, shell=True)

if __name__ == "__main__":
	path_to_aug = sys.argv[1]
	path_to_save = sys.argv[2]
	species = sys.argv[3]
	genbankF = sys.argv[4]
	percentage = int(sys.argv[5])

	#Creation of new directory
	speciesfolder = newDirectory(path_to_aug, species)

	#Copy and edit all files from 'standard' predefined files
	copyGenericFiles(path_to_aug, species, speciesfolder)
	changeFilenames(speciesfolder,species)
	changeParameterFile(speciesfolder, species)

	#Creation of training and test set of genbankfile
	train, test = splitGenbankInTrainTest(genbankF, percentage)
	writeGB(path_to_save, train, test)

	#Optimise augustus parameters
	optimiseParameters(path_to_aug, species, path_to_save)
	etraining(path_to_aug, species, path_to_save)
	
	#Check accuracy of your prediction
	
