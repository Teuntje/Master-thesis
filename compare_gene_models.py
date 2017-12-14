#!/usr/bin/env python

"""
Author: Teuntje Peeters
Student number: 920301645120

Script for: comparing a reference transcript to a predicted transcript 
by using the tool augustus

This script only compares positions and it calculates the recall and precision
It was written to get an idea of comparing gene models

"""

from sys import argv
import subprocess
import os.path

def runAugustus(fasta, genemodels="complete", species="saccharomyces_cerevisiae_S288C", output="output_augustus.gff"):
	"""Run augustus on fastafile

	fasta - fasta file with genome
	"""
	if os.path.isfile(output):
		pass
	else:
		cmd = "augustus --genemodel={} --species={} {} > {}".format(genemodels, species, fasta, output)
		print cmd
		e = subprocess.check_call(cmd, shell=True)
	return

def parse_referenceGFF(gff):
	"""Parse gff file
	
	gff - gff file
	"""
	mRNA = []
	with open(gff) as inFile: 
		for line in inFile:
			if line.startswith("#"):
				pass
			else:
				if "mRNA" in line:
					start = line.split()[3]
					stop = line.split()[4]
					p = line.split()[8].split(";")
					for i in p:
						if "Parent" in i:
							parent = i.split("=")[1]
							if "A" in parent:
								pass
							else:
								mRNA.append([start, stop, parent])
	return  mRNA

def parse_augustusGFF(output):
	"""Parse augustus output
	
	output from augustus as gff file
	"""
	transcript = []
	with open(output) as inFile:
		for line in inFile:
			if line.startswith("#"):
				pass
			else:
				if "transcript" == line.split()[2]:
					start = line.split()[3]
					stop = line.split()[4]
					ID = line.split()[8]
					transcript.append([start, stop, ID])
	return transcript

def compareGeneModels(reference, augustus):
	"""Compare gene models

	reference - reference annotation
	augustus - prediction of gene models by augustus
	"""
	tp = []
	fp = []
	for r in reference:
		for a in augustus:
			if r[0] == a[0] and r[1] == a[1]:
				tp.append([r[0], r[1], a[2], r[2]])

	#get false negatives
	genes_tp = [a[3] for a in tp]
	genes_ref = [r[2] for r in reference]

	false_negatives = []
	for i in genes_ref:
		if i not in genes_tp:
			false_negatives.append(i)
	#print false_negatives

	#get false positives
	false_positives=[]
	ID_tp = [a[2] for a in tp]
	ID_aug = [a[2] for a in augustus]
	for i in ID_aug:
		if i not in ID_tp:
			false_positives.append(i)

	return tp, false_negatives, false_positives

def calculation(tp, fn, fp):
	"""Calculaten the recall and precision

	tp - true positives
	fn - false negatives
	fp - false positives
	"""
	recall = float(len(tp))/float((len(tp)+len(fn)))
	precision = float(len(tp))/float((len(tp)+len(fp)))
	return recall, precision

def outputfile(tp):
	"""create outputfile

	tp - true positive output
	"""
	output = open("outputfile.txt", "w")
	output.write("Start \t Stop \t AugID \t GeneID")
	for i in tp:
		output.write("\n")
		output.write("{} \t {} \t {} \t {}".format(i[0], i[1], i[2], i[3]))

if __name__ == "__main__":
	fasta = argv[1]
	gff = argv[2]
	#runAugustus(fasta, genemodels="complete", species="saccharomyces_cerevisiae_S288C", output="output_augustus.gff")
	referenceGFF = parse_referenceGFF(gff)
	augustusGFF = parse_augustusGFF(output="output_augustus.gff")
	#tp, fn, fp = compareGeneModels(referenceGFF, augustusGFF)
	#recall, precision = calculation(tp, fn, fp)
	print "Recall: ", round(recall, 2)
	print "Precision: ", round(precision, 2)
	outputfile(tp)
