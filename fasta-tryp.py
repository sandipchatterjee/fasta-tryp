#! /usr/bin/python

import sys

if len(sys.argv) != 2:
	print 'Need one argument -- a filename'
	print 'Proper usage: >>fastaparse.py somefile.fasta'
	sys.exit()

fileName = sys.argv[1]
fileDict = {}
inputFile = open(fileName)

for fileLine in inputFile:
	if '>' in inputFile.readline():
		

inputFile.close()