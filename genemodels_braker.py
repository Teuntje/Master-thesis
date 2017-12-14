#!/usr/bin/env python

"""
Author: Teuntje Peeters
Student number: 920301645120

Script for: comparing a reference transcript to a predicted transcript 
by using the tool augustus
This script is not used in further analysis. 
The purpose was to get familiar with comparing genemodels and getting some
output

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

def parseGFF2(gff):
	"""Parse gff file

	Obtain the transcripts, introns and exons

	gff - gff file
	"""
	#Initiate three different dictionaries
	genes = {}
	introns = {}
	exons = {}

	#Open and read the file
	with open(gff) as inFile:
		for line in inFile:
			#If line does not start with #
			if not line.startswith("#"):
				#If transcript in line, save start and stop position
				#Save the key for the dictionary (gene ID)
				#Save the percentage and add to the dictionary
				if "transcript" == line.split()[2]:
					start = line.split()[3]
                                        stop = line.split()[4]
                                        key = line.split()[8]
					percentage = line.split()[5]
					genes[key] = [[start, stop, percentage]]
					introns[key] = []
					exons[key] = []
				#If intron in line save the start, stop position and percentage
				#Save the key for the dictionary
				if "intron" == line.split()[2]:
					i_start = line.split()[3]
					i_stop = line.split()[4]
					percentage = line.split()[5]
					key = line.split()[9].split(";")[0][1:-1]
					#Check if key does not start with a number or add values to the dictionary
					if key[0].isdigit():
						pass
					else:
						introns[key].append([i_start, i_stop, percentage])
				#If CDS is in line, save start, stop and percentage
				#Save the key for the dictionary
				if "CDS" == line.split()[2]:
					e_start = line.split()[3]
					e_stop = line.split()[4]
					percentage = line.split()[5]
					key = line.split()[9].split(";")[0][1:-1]
					if not key[0].isdigit():
						exons[key].append([e_start, e_stop, percentage])

	#Check for empty values and delete these
	genes = emptyValues(genes)
	exons = emptyValues(exons)
	introns = emptyValues(introns)

	return genes, exons, introns

def emptyValues(dic):
	"""Check and remove empty values in dictionary

	dic - dictionary
	"""
	for key in list(dic.keys()):
	    if dic[key] == []:
        	del dic[key]
	return dic


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

def compareGeneModels2(braker, augustus):
	"""Compare gene models to return start, stop and complete matches
	
	- braker can be dictionary with transcripts, introns and exons
	- augustus can be dictionary with transcripts, introns and exons
	"""
	#Initiate lists for matches, start position matches and stop position matches
	matches = []
	start_match = []
	stop_match = []
	#Compare everything with everything
	for bk, bv in braker.iteritems():
		for ak, av, in augustus.iteritems():
			if len(bv) > 1 or len(av) > 1:
				for b in bv:
					for a in av:
						if b[0] == a[0] and b[1] == a[1]:
							matches.append([b[0], b[1], bk, ak])
						elif b[0] == a[0] and b[1] != a[1]:
							start_match.append([b[0], b[1], a[1], bk, ak])
						elif b[0] != a[0] and b[1] == a[1]:
							stop_match.append([b[0], b[1], a[1], bk, ak])
			else:
				if bv[0][0] == av[0][0] and bv[0][1] == av[0][1]:
					matches.append([bv[0][0], bv[0][1], bk, ak])
				elif bv[0][0] == av[0][0] and bv[0][1] != av[0][1]:
					start_match.append([bv[0][0], bv[0][1], av[0][1], bk, ak])
				elif bv[0][0] != av[0][0] and bv[0][1] == av[0][1]:
					stop_match.append([bv[0][0], bv[0][1], av[0][0], bk, ak])
	return matches, start_match, stop_match

def calculation(tp, fn, fp):
	"""Calculaten the recall and precision

	tp - true positives
	fn - false negatives
	fp - false positives
	"""
	recall = float(len(tp))/float((len(tp)+len(fn)))
	precision = float(len(tp))/float((len(tp)+len(fp)))
	return recall, precision

def outputfile(tp, filename):
	"""create outputfile

	tp - output to be written in a file
	filename - string with filename
	"""
	output = open("{}.txt".format(filename), "w")
	output.write("Start \t \t Stop \t \t  extra \t \t brakerID \t \t augustusID")
	for i in tp:
		if len(i) > 4:
			output.write("\n")
			output.write("{} \t \t {} \t \t {} \t \t {} \t \t {}".format(i[0], i[1], i[2], i[3], i[4]))
		else:
			output.write("\n")
			output.write("{} \t \t {} \t \t {} \t \t {}".format(i[0], i[1], i[2], i[3]))

if __name__ == "__main__":
	braker = argv[1]
	augustus = argv[2]
	#runAugustus(fasta, genemodels="complete", species="saccharomyces_cerevisiae_S288C", output="output_augustus.gff")
	brakerTranscripts, brakerExons, brakerIntrons = parseGFF2(braker)
	print "Amount of braker transcripts:", len(brakerTranscripts)
	print "Amount of braker exons:", len(brakerExons)
	print "Amount of braker introns:", len(brakerIntrons)
	augustusTranscripts, augustusExons, augustusIntrons = parseGFF2(augustus)
	print "Amount of augustus transcripts:", len(augustusTranscripts)
	print "Amount of augustus exons:", len(augustusExons)
	print "Amount of augustus introns:", len(augustusIntrons)
	print "start transcript comparison"
	transcriptMatches, transcriptStartMatches, transcriptStopMatches = compareGeneModels2(brakerTranscripts, augustusTranscripts)
	print "Length transcript matches:", len(transcriptMatches)
	outputfile(transcriptMatches, "transcript_matches"); outputfile(transcriptStartMatches, "transcript_start_matches"); outputfile(transcriptStopMatches, "transcript_stop_matches")
	print "transcript comparison, done"
	print "Outputfiles made"
	print "Continue to exon comparison .."
	exonMatches, exonStartMatches, exonStopMatches = compareGeneModels2(brakerExons, augustusExons)
	outputfile(exonMatches, "exon_matches"); outputfile(exonStartMatches, "exon_start_matches"); outputfile(exonStopMatches, "exon_stop_matches")
	print "exon comparison, done"
	print "Outputfiles made"
	print "continue to intron comparison"
	intronMatches, intronStartMatches, intronStopMatches = compareGeneModels2(brakerIntrons, augustusIntrons)
	outputfile(intronMatches, "intron_matches"); outputfile(intronStartMatches, "intron_start_matches"); outputfile(intronStopMatches, "intron_stop_matches")
	print "Outputfiles made" 
	print "intron comparison, done"
#	recall, precision = calculation(tp, fn, fp)
#	print "Recall: ", round(recall, 2)
#	print "Precision: ", round(precision, 2)
#	outputfile(transcriptStopMatches, "transcript_stop_matches")
