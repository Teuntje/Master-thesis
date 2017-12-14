#!/usr/bin/env python

"""
Author: Teuntje Peeters
Student number: 920301645120

Script for: selection of the best gene models to create the final
gene model file

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
			key = line.split()[0]
			if key in genes:
				genes[key].append(line.split())
			else:
				if float(line.split()[2]) >= 70.0:
					genes[key] = [line.split()]
	return genes

def gffParser_without_matches(gff, bbh, filename, o_gm):
	"""Parse gff file
	"""
	exportfile = open("{}_without_matches.gff".format(filename), "w")
	bbh = bbh.keys()
	with open(gff) as inFile:
		for line in inFile:
			if not line.startswith("#"):
				for hit in bbh:
					if hit in line:
						exportfile.write(line)
					elif hit.split(".")[0] == line.split()[8] and line.split()[2] == "gene":
						exportfile.write(line)
				for nm in o_gm:
					nm1 = "{}.t1".format(nm)
					if nm in line and line.split("\t")[2] == "gene":
						exportfile.write(line)
					if nm1 in line:
						exportfile.write(line)

def gffParser_only_matches(gff, bbh, filename):
        """Parse gff file
        """
        exportfile = open("{}_only_matches.gff".format(filename), "w")
        bbh = bbh.keys()
        with open(gff) as inFile:
                for line in inFile:
                        if not line.startswith("#"):
                                for hit in bbh:
                                        if hit in line:
                                                exportfile.write(line)
                                        elif hit.split(".")[0] == line.split()[8] and line.split()[2] == "gene":
                                                exportfile.write(line)

def extractGenemodels(gff):
	"""Extract gene models from original gff file
	"""
	genemodels = []
	with open(gff) as inFile:
		for line in inFile:
			if not line.startswith("#"):
				if line.split("\t")[2] == "transcript":
					line = line.strip().split("\t")
					if line[8].endswith(".t1"):
						genemodels.append(line[8].split(".")[0])
	return genemodels

def getGenemodelsWithoutMatches(bbh, original_genemodels):
	"""Get the genemodels without matches

	"""
	bbh1 = []
	for i in bbh:
		bbh1.append(i.split(".")[0])

	return list(set(original_genemodels) - set(bbh1))

if __name__ == "__main__":
	#arvg1 = BBH query vs subject
	query_vs_sub = argv[1]
	#argv2 = genemodels of augustus
	gff = argv[2]
	original_genemodels = extractGenemodels(gff)
	#argv3 = filename of output file
	bbh = csvparser(query_vs_sub)
	no_match_gm = getGenemodelsWithoutMatches(bbh, original_genemodels)
	gffParser_only_matches(gff, bbh, argv[3])
	gffParser_without_matches(gff, bbh, argv[3], no_match_gm)
