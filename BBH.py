#!/usr/bin/env python

"""
Author: Teuntje Peeters
Student number: 920301645120

Script for selection of the Best bidirectional hits and adding
the lengths of the sequences + calculates the Marnix index

Arguments:
arg1 = organism1_vs_organism2 = blast output in csv format
arg2 = organism2_vs_organism1 = blast output in csv format
arg3 = organism1 fasta sequences
arg4 = organism2 fasta sequences 

"""

from sys import argv
import subprocess
import os.path
from textwrap import wrap
import time

def csvparser(csvfile):
	"""Parse csv

	"""
	with open(csvfile) as inFile:
		genes = {}
		for line in inFile:
			key = line.split()[0].split(".")[0]
			if key in genes:
				genes[key].append(line.split())
			else:
				genes[key] = [line.split()]
	return genes

def getAmountofHits(dict):
	"""Get the amount of hits in total

	"""
	count = 0
	for key, value in dict.items():
		for v in value:
			count += 1
	return count

def getBestHits(csv_dict):
	"""Filter best hits

	"""
#	print csv_dict
	bestHits = {}
	for k, v in csv_dict.items():
		maxBIT = max([float(val[11]) for val in v])
		bestHits[k] = [val for val in v if str(maxBIT) == str(float(val[11]))][0]
	return bestHits

def getBBH(dict1, dict2):
	"""Get the best bidirectional hit

	"""
	bbh = []
	bbh_terreus_vs_oryzae_braker = []
	bbh_oryzae_vs_terreus_augustus = []
	for k, v in dict1.items():
		for key, value in dict2.items():
			if [k, v[1].split(".")[0]] == [value[1].split(".")[0], key]:
				#print k, v[1].split(".")[0], "\t", key, value[1].split(".")[0]
				bbh.append(v)
				bbh.append(value)
				bbh_oryzae_vs_terreus_augustus.append(v)
				bbh_terreus_vs_oryzae_braker.append(value)
	return bbh, bbh_oryzae_vs_terreus_augustus, bbh_terreus_vs_oryzae_braker

def writefile(BBH, filename):
	"""Write csv file

	"""
	if os.path.isfile("{}.csv".format(filename)):
		pass
	else: 
		with open("{}.csv".format(filename), "w") as inFile:
			for line in BBH:
				for pos in line:
					inFile.write("{}\t".format(pos))
				inFile.write("\n")

def fastaParser(fastafile):
	"""Get for each genemodel the length of the sequence

	"""
	fastaDict = {}
	with open(fastafile) as inFile:
		for line in inFile:
			if line.startswith(">"):
				key = line.split(">")[1].strip()
				fastaDict[key] = []
			else:
				fastaDict[key].append(line.strip())
	fastaDict = {key:len("".join(value)) for key, value in fastaDict.items()}
	return fastaDict

def addLengthSequence(csv, fasta1, fasta2):
	"""Add length to sequences in csv dictionary

	"""
	for k, v in csv.items():
		csv[k] = [item for sublist in v for item in sublist]
	for k, v in csv.items():
		csv[k].append(fasta1[v[0]])
		csv[k].append(fasta2[v[1]])
	return csv.values()

def calcMarnixIndex(csv):
	"""Calculate a variant of the jaccard index

	"""
	new_csv = []
	for line in csv:
		marnix_index= ((float(line[9])-float(line[8]))+(float(line[7])-float(line[6])))/(float(line[12])+float(line[13]))
		line.append(repr(round(marnix_index, 3)))
		new_csv.append(line)
	return  new_csv

if __name__ == "__main__":
	start_time = time.time()

	#Get right datasets
	ory_vs_ter = csvparser(argv[1])
	ter_vs_ory = csvparser(argv[2])
	ory_fasta = argv[3]
	ter_fasta = argv[4]
	#Get the gene length in dictionary
	ory_gene_length = fastaParser(ory_fasta)
	ter_gene_length = fastaParser(ter_fasta)
	#Some statistics
	print "Amount of Blast hits Oryzae vs Nidulans:", getAmountofHits(ory_vs_ter)
	print "Amount of Blast hits Nidulans vs Oryzae:", getAmountofHits(ter_vs_ory)
	#Get the best hits for both files

	bestOryVSTer = getBestHits(ory_vs_ter)
	bestTerVSOry = getBestHits(ter_vs_ory)
	#print bestOryVSTer
	#Again, some statistics
	print "Amount of best hits Oryzae vs Nidulans:", len(bestOryVSTer)
	print "Amount of best hits Nidulans vs Oryzae:", len(bestTerVSOry)
	#Get the best bidirectional hits
	BBH, BBH_oryzae_vs_terreus_augustus, BBH_terreus_vs_oryzae_braker = getBBH(bestOryVSTer, bestTerVSOry)
	#Again some statistics
	print "Amount of BBH:", len(BBH)
	print "Amount of BBH pairs:", len(BBH)/2

	writefile(BBH_oryzae_vs_terreus_augustus, "BBH_oryzae_vs_nidulans_augustus")
	writefile(BBH_terreus_vs_oryzae_braker, "BBH_nidulans_vs_oryzae_braker")
	#Run this one only if you don't want to run the above (to avoid time)
	BBH_oryzae_vs_terreus_augustus = csvparser("BBH_oryzae_vs_nidulans_augustus.csv")
	BBH_terreus_vs_oryzae_braker = csvparser("BBH_nidulans_vs_oryzae_braker.csv")
	
	BBH_oryzae_vs_terreus_augustus = addLengthSequence(BBH_oryzae_vs_terreus_augustus, ory_gene_length, ter_gene_length)
	BBH_terreus_vs_oryzae_braker = addLengthSequence(BBH_terreus_vs_oryzae_braker, ter_gene_length, ory_gene_length)
	
	BBH_oryzae_vs_terreus_augustus = calcMarnixIndex(BBH_oryzae_vs_terreus_augustus)
	BBH_terreus_vs_oryzae_braker = calcMarnixIndex(BBH_terreus_vs_oryzae_braker)
	
	writefile(BBH_oryzae_vs_terreus_augustus, "BBH_oryzae_vs_nidulans_augustus_lengths")
	
	writefile(BBH_terreus_vs_oryzae_braker, "BBH_nidulans_vs_oryzae_braker_lengths")
	print "Time it took: ", time.time() - start_time, "to run"
