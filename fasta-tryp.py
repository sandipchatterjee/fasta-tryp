#! /usr/bin/python

##	fasta-tryp.py
##	For generating a nonredundant peptide database in FASTA format from an input FASTA-formatted protein database
##	Sandip Chatterjee

import sys
import pprint	##	Pretty Print

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
############# need to account for presence of proline #############

##	Info from Expasy:
##	Preferentially cleaves at Arg and Lys in position P1 with higher rates for Arg (Keil, 1992), especially at high pH (but treated equally in the program). Pro usually blocks the action when found in position P1', but not when Lys is in position P1 and Trp is in position P2 at the same time. This blocking of cleavage exerted by Pro in position P1' is also negligible when Arg is in position P1 and Met is in position P2 at the same time (other reports say that the block exhibited by Pro can be circumvented by Glu being in P2). 

############# can probably filter out peptides below a certain length (maybe... 6 amino acids?) #############
def trypsinDigest(protSeq):
	
	trypticPeptides = protSeq.replace('K','K_').replace('R','R_').split('_')

	trypticPeptides = filter(lambda x: len(x)>6,trypticPeptides)			##	only keep tryptic peptides with length > 6

	return trypticPeptides

###################################################################################################

##	sample protein data
sampleProtein = "MSGSAAETQAPEQISGGWTPVVEETPAPPAPPVRRRGVPAGPVIIAAGEAGVMAGSSLYAAAGVPGLIAAGAVAGTAAAARTVRVAQRRAASRASTSGGGAGRGRSGSGPVPAALRRLAGAGRRRGAGGGPGTGKGAGTGTGTGRRRRGAGTGRRPGSGRAGGAAGRTAGGRRAGGGPLTRAARGTGGQAARAARRSGAAAVRAARLRGAGMVRAARIGRSAAGRLAAAVRQARSPWAAARAARRGRYRMVAGQRRRRGGGLGRLVRTGLGLLAALAAVGWWLGHRGWRLTRAGWRRLRGPRPGPDGRPGTDLIGPLAQLRTPPIGGPAAPVPPPAGGGAMDPQAILAGARAALAGVAFATFGRLPEMPEMPPIPGGPVASLMPGHDRGDAAMAARDRDGRSLAGGSGGVFLLPEQAETMMQTAAAYAPPSGPQIGRDMAQVPQAVEHIAAAIAGLSRLAADELPVHETIRELLDSAAAQIVAAAGTLANLGPLFELVHQTDLHRLRQPRPGEGLWDSHHAGDA"
samplePeptides = trypsinDigest(sampleProtein)
sampleFastaPeptides = [">gi|13449237|ref|NP_085453.1| hypothetical protein pFQ12_p05 [Frankia sp. CpI1]", samplePeptides]

sampleFastaRecord = [sampleFastaPeptides[0],sampleProtein]

###################################################################################################

##	function parseFastaRecord(fastaRecord)
##	input: a single FASTA record from a FASTA protein file (two-item list)
##	output: a two-item list with the format [{'NCBIID':VALUE, 'proteinName':VALUE, 'organismID':VALUE, 'protSeq':VALUE}, [PEPTIDE1, PEPTIDE2, PEPTIDE3, ...]]
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

	parsedFastaRecord = [{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protSeq': protSeq}, trypsinDigest(protSeq)]

	return parsedFastaRecord

###################################################################################################

##	function registerPeptides
##	Registers peptides from a single FASTA record in allPeptides dictionary
##	input: FASTA record (modified; Python list)
##	input list format: [{'NCBIID':VALUE, 'proteinName':VALUE, 'organismID':VALUE, 'protSeq':VALUE}, [PEPTIDE1, PEPTIDE2, PEPTIDE3, ...]]
##	output: allPeptides dictionary with added peptide entries
##	output list format: {'PEPTIDESEQUENCE':[{'NCBIID': NCBIID, 'proteinName': proteinName, 'organismID': organismID, 'protSeq': protSeq}, ...]}
##	output list format: Dictionary in which each Peptide sequence (key) is mapped to a list of length >= 1
##	output list format: Each list item is a dictionary of information containing the source of the peptide
##	output list format: One peptide may map to one or more proteins
############# need to see what data should be extracted from the first > line of the FASTA entry... #############
############# eventually need to have more information about each key (more than just parent protein) -- start site in protein sequence, bacterial species/organism ID, etc. #############
def registerPeptides(fastaRecord):
	for peptide in fastaRecord[1]:
		if peptide not in allPeptides:
				allPeptides[peptide] = [fastaRecord[0]]		##	creating a list of dictionaries for each peptide in allPeptides dictionary
		else:											## in this case, need to *append* the PROTEIN INFO DICTIONARY to allPeptides[peptide]
			print "Redundant peptide"
#			print fastaRecord[0]		#remove

	print "END"

	return allPeptides

###################################################################################################


##	function readFasta()
##	Reads in all entries from a protein FASTA-format file
##	Calls function parseFastaRecord() for each imported record
##	Calls registerPeptides() for each parsedFastaRecord
##	input: FASTA-formatted record file
##	output: complete allPeptides dictionary
def readFasta():
	fileName = sys.argv[1]
	fileList = []
	fastaRecord = []
	inputFile = open(fileName, "r")

	for fileLine in inputFile:					##	read entire file into list... probably a bad idea
		fileList.append(fileLine)

	for fileLine in fileList:
		if '>' in fileLine:						##	store each FASTA Record in 2-item list "fastaRecord" temporarily
			fastaRecord.append(fileLine)
	
			print "fastaRecord:"
			print fastaRecord



#		parsedFastaRecord = parseFastaRecord(fastaRecord)	##	then call parseFastaRecord(fastaRecord)

	inputFile.close()

	return allPeptides

###################################################################################################

#parsedFastaRecord = parseFastaRecord(sampleFastaRecord)		## remove later
#allPeptides = registerPeptides(parsedFastaRecord)

###################################################################################################

readFasta()

print "----------------------------------------------------------"
#print "allPeptides dictionary: "
#print allPeptides

pp = pprint.PrettyPrinter()
pp.pprint(allPeptides)

print "length of dictionary: "
print len(allPeptides)

print "----------------------------------------------------------"
print allPeptides

#registerPeptides(sampleFastaPeptides)

###################################################################################################

