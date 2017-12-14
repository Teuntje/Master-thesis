#!/usr/bin/env python

"""
Author: Teuntje Peeters
Student number: 920301645120

Script for: converting gff to fasta file

"""

from sys import argv
import subprocess
import os.path
from textwrap import wrap

def parseGff(gff):
	"""Parse gff
	"""
	output = {}
	keepreading = False
	with open(gff) as inFile: 
		for line in inFile:
			if not line.startswith("#"):
				keepreading = False
				if "transcript" == line.split()[2]:
					key = line.split()[8]
					output[key] = []
			else:
				if "protein sequence" in line:
					if "]" in line:
						output[key].append(line.split()[4][1:].split("]")[0])
					else:
						output[key].append(line.split()[4][1:])
				elif line.isupper():
					if "]" in line:
						line = line.split("]")[0]
						output[key].append(line.split()[1])
					if ":" in line:
						continue
					else:
						output[key].append(line.split()[1])
	output = {key: "".join(value) for key, value in output.items()}
	return output

def writeFasta(dict, filename):
	"""Write fasta file
	"""
	output = open("{}".format(filename), "w")
	for key, value in dict.items():
		output.write(">{}\n".format(key))
		output.write('\n'.join(['\n'.join(wrap(line, width=45)) for line in value.splitlines()]))
		output.write('\n')

if __name__ == "__main__":
	gff = argv[1]
	filename = argv[2]
	gff_dict = parseGff(gff)
	writeFasta(gff_dict, filename)
