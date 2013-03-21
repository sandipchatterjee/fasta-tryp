#! /usr/bin/python

##	fasta-tryp.py
##	Sandip Chatterjee
##
##	Usage:
##		>>python fasta-tryp.py proteome_fasta_file.fasta
##
##	For generating a nonredundant peptide database in FASTA format from an input FASTA-formatted protein database
##	1) reads in all proteins from input file (standard FASTA format)
##	2) generates tryptic peptides from each protein based on PeptideCutter rules
##	3) generates nonredundant dictionary of tryptic peptides of length > 6
##	4) outputs a nonredundant peptide file for use in spectral matching searches
##	5) outputs an index file with parent proteins and organisms for each peptide 
##	(all output on one line; can be searched using grep '>peptide1234|' < *.fastaindex)

######	to do:
######	account for half-tryptic peptides?

import sys

if len(sys.argv) != 2:
	print 'Need one argument -- a filename'
	print 'Proper usage: >>fasta-tryp.py somefile.fasta'
	sys.exit()

###################################################################################################

#global peptide dictionary
allPeptides = {}

###################################################################################################

##	function trypsinDigest
##	input: protein sequence (string, case-insensitive)
##	output: tryptic peptides from input protein sequence (list of peptides, upper-case)
##	
##	Info from Expasy Peptide Cutter (http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#Tryps):
##	Preferentially cleaves at Arg and Lys in position P1 with higher rates for Arg (Keil, 1992), especially at high pH (but treated equally in the program). 
##	
##	SPECIAL CASES FOR PROLINE AT P1':
##	Pro usually blocks the action when found in position P1', 
##	but not when Lys is in position P1 and Trp is in position P2 at the same time. 
##	This blocking of cleavage exerted by Pro in position P1' is also negligible when Arg is in position P1 and Met is in position P2 at the same time 
##	
##	SPECIAL CASES FOR LYSINE AT P1:
##	Furthermore, if Lys is found in position P1 the following situation considerably block the action of trypsin: 
##	Either Asp in position P2 and Asp in position P1' or 
##	Cys in position P2 and Asp in position P1' or 
##	Cys in position P2 and His in position P1' or 
##	Cys in position P2 and Tyr in position P1'. 
##	
##	SPECIAL CASES FOR ARGININE AT P1:
##	A likewise considerable block of trypsin action is seen when Arg is in P1 and the following situations are found: 
##	Either Arg in position P2 and His in position P1' or 
##	Cys in position P2 and Lys in position P1' or 
##	Arg in position P2 and Arg in position P1'.


def trypsinDigest(protSeq):
	
	trypticProtSeq = protSeq.upper().replace('K','K_').replace('R','R_')	##	define all Lys, Arg potential cleavage sites as '_'

	##	special cases for trypsin cleavage with proline at P1'
	trypticProtSeq = trypticProtSeq.replace('_P','P')		##	if P1' residue (C-terminal to cleavage site) is Proline, remove cleavage site
	trypticProtSeq = trypticProtSeq.replace('WKP','WK_P')	##	if P1 is K and P2 is W, Proline at P1' will not block cleavage
	trypticProtSeq = trypticProtSeq.replace('MRP','MR_P')	##	if P1 is R and P2 is M, Proline at P1' will not block cleavage

	##	special cases for trypsin cleavage with lysine at P1
	trypticProtSeq = trypticProtSeq.replace('DK_D','DKD')		##	if P1 is K, P2 is D, and P1' is D, remove cleavage site
	trypticProtSeq = trypticProtSeq.replace('CK_D','CKD')		##	if P1 is K, P2 is C, and P1' is D, remove cleavage site
	trypticProtSeq = trypticProtSeq.replace('CK_H','CKH')		##	if P1 is K, P2 is C, and P1' is H, remove cleavage site
	trypticProtSeq = trypticProtSeq.replace('CK_Y','CKY')		##	if P1 is K, P2 is C, and P1' is Y, remove cleavage site

	##	special cases for trypsin cleavage with arginine at P1
	trypticProtSeq = trypticProtSeq.replace('R_R_H','R_RH')		##	if P1 is R, P2 is R, and P1' is H, remove cleavage site
	##	want to keep R_RH but not RR_H
	trypticProtSeq = trypticProtSeq.replace('CR_K','CRK')		##	if P1 is R, P2 is C, and P1' is K, remove cleavage site
	trypticProtSeq = trypticProtSeq.replace('R_R_R','R_RR')		##	if P1 is R, P2 is R, and P1' is R, remove cleavage site
	##	want to keep R_RR but not RR_R

	trypticPeptides = trypticProtSeq.split('_')		##	cleave at all sites in protSeq marked by '_'

	trypticPeptides = filter(lambda x: len(x)>6,trypticPeptides)			##	only keep tryptic peptides with length > 6

	return trypticPeptides

###################################################################################################

##	function parseFastaRecord(fastaRecord)
##	input: a single FASTA record from a FASTA protein file (two-item list)
##	output: a two-item list with the format [{'NCBIID':VALUE, 'proteinName':VALUE, 'organismID':VALUE, 'protSeq':VALUE, 'protPosition': -2}, [PEPTIDE1, PEPTIDE2, PEPTIDE3, ...]]
##	output details: calls trypsinDigest() function using given protein sequence from fastaRecord

#################### This probably isn't the best way to split the informational line... using Split('[') -- because not every entry has brackets ######################
def parseFastaRecord(fastaRecord):
	
	parsedFastaRecord = [{},[]]
	infoLine = fastaRecord[0]

	##	parsing infoLine...

	##	...for NCBI Protein ID
	splitInfoLine1 = infoLine.split('|')
	NCBIID = splitInfoLine1[3]

	##	...for protein name
	splitInfoLine1 = infoLine.split('|')
	proteinName = splitInfoLine1[4].split('[')[0]	##	remove Organism info (if present)

	##	...for organism ID
	splitInfoLine1 = infoLine.split('[')
	if ']' in splitInfoLine1[len(splitInfoLine1)-1]:	
		organismID = splitInfoLine1[len(splitInfoLine1)-1].strip(']')	## working for now, but need to account for entries without brackets in infoLine []
	else:
		organismID = None

	##	Obtaining protein sequence
	protSeq = fastaRecord[1]

	## construct dictionary

	parsedFastaRecord = [{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protSeq': protSeq, 'protPosition': -2}, trypsinDigest(protSeq)]

	return parsedFastaRecord

###################################################################################################

##	function registerPeptides
##	Registers peptides from a single FASTA record in allPeptides dictionary
##	input: FASTA record (modified; Python list)
##	input list format: [{'NCBIID':VALUE, 'proteinName':VALUE, 'organismID':VALUE, 'protSeq':VALUE}, [PEPTIDE1, PEPTIDE2, PEPTIDE3, ...]]
##	output: allPeptides dictionary with added peptide entries
##	output list format: {'PEPTIDESEQUENCE':[{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protPosition': INT}, ...]}
##	output list format: Dictionary in which each Peptide sequence (key) is mapped to a list of length >= 1
##	output list format: Each list item is a dictionary of information containing the source of the peptide
##	output list format: One peptide may map to one or more proteins
############# protPosition doesn't work properly when there are two identical tryptic peptides in a protein sequence #############
############# but no way to distinguish these by mass spec anyway... 												 #############
def registerPeptides(fastaRecord):

	for peptide in fastaRecord[1]:
		if peptide not in allPeptides:
			allPeptides[peptide] = []		##	initialize list of protein records for new peptide

		NCBIIDTemp = fastaRecord[0]['NCBIID']
		proteinNameTemp = fastaRecord[0]['proteinName']
		organismIDTemp = fastaRecord[0]['organismID']
		protPositionTemp = fastaRecord[0]['protSeq'].find(peptide)

		allPeptides[peptide].append({'NCBIID': NCBIIDTemp,'proteinName': proteinNameTemp, 'organismID': organismIDTemp,'protPosition': protPositionTemp})

	return allPeptides

###################################################################################################


##	function readFasta()
##	Reads in all entries from a protein FASTA-format file (one entry at a time)
##	Calls function parseFastaRecord() for each imported record
##	Calls registerPeptides() for each parsedFastaRecord
##	input: FASTA-formatted record file
##	output: complete allPeptides dictionary
def readFasta():
	fileName = sys.argv[1]
	fileList = []
	fastaRecord = []
	inputFile = open(fileName, "r")
	
	outputFile = open('reformatted_'+fileName, "w")

	for line in inputFile:				##	generates a properly formatted FASTA file. InfoLine ('>') on one line, entire sequence on next line
		if '>' in line:
			outputFile.write('\n'+line)		##	keep in mind that this prints a blank line at the beginning of the reformatted file
		else:
			outputFile.write(line.strip('\n').strip('\t'))


	outputFile.close()
	inputFile.close()

	inputFile = open('reformatted_'+fileName, "r")	##	read in properly formatted FASTA file (one line per protein sequence)
	inputFile.readline()							##	read blank line at beginning of reformatted file

	while True:
		line1 = inputFile.readline().strip('\n')
		line2 = inputFile.readline().strip('\n')
		if not line1:
			break
		fastaRecord.append(line1)	##	info line of FASTA record ('>')
		fastaRecord.append(line2)	##	sequence line of FASTA record
		allPeptides = registerPeptides(parseFastaRecord(fastaRecord))
		fastaRecord = []

	inputFile.close()

	return allPeptides

###################################################################################################

##	function outputFasta(allPeptides)
##	Outputs contents of allPeptides dictionary to FASTA formatted-file
##	input: complete allPeptides dictionary
##	output: FASTA file and index file

##	format of FASTA file:
##	>PEPTIDE1
##	TRYPTICPEPTIDESEQUENCE1
##	>PEPTIDE2
##	TRYPTICPEPTIDESEQUENCE2
##	...

##	format of index file:
##	>PEPTIDE1
##	[{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protPosition': INT}]
##	>PEPTIDE2
##	[{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protPosition': INT},
##	{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protPosition': INT}]
##	>PEPTIDE3
##	[{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protPosition': INT}]
##	>PEPTIDE4
##	[{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protPosition': INT},
##	{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protPosition': INT},
##	{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protPosition': INT}]
##	...

def outputFasta(allPeptides):

	##	generate FASTA file
	outputFastaFile = open(sys.argv[1].replace('.fasta','')+'_peptides_nonredundant'+'.fasta', 'w')

	peptideCount = 1
	for key in allPeptides:
		outputFastaFile.write('>'+'peptide'+str(peptideCount)+'\n'+key+'\n')
		peptideCount += 1

	print "Generated file "+sys.argv[1].replace('.fasta','')+'_peptides_nonredundant'+'.fasta'
	outputFastaFile.close()

	##	generate index file
	outputIndexFile = open(sys.argv[1].replace('.fasta','')+'_peptides_nonredundant'+'.fastaindex', 'w')

	peptideCount = 1
	for key in allPeptides:
		outputIndexFile.write('>'+'peptide'+str(peptideCount)+'|')	## pipe character after peptide ID
		peptideCount += 1
		for item in range(len(allPeptides[key])):
			if allPeptides[key][item]['organismID']:
				outputIndexFile.write(key+','+allPeptides[key][item]['NCBIID']+','+allPeptides[key][item]['proteinName']+','+allPeptides[key][item]['organismID']+','+str(allPeptides[key][item]['protPosition'])+'|')
			else:
				outputIndexFile.write(key+','+allPeptides[key][item]['NCBIID']+','+allPeptides[key][item]['proteinName']+','+''+','+str(allPeptides[key][item]['protPosition'])+'|')
		outputIndexFile.write('\n')

	print "Generated file "+sys.argv[1].replace('.fasta','')+'_peptides_nonredundant'+'.fastaindex'
	outputIndexFile.close()

###################################################################################################

readFasta()
outputFasta(allPeptides)

print "Finished"

###################################################################################################