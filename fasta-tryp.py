#! /usr/bin/python

##	fasta-tryp.py
##	For generating a nonredundant peptide database in FASTA format from an input FASTA-formatted protein database
##	Sandip Chatterjee

###################################################################################################

#global peptide dictionary
allPeptides = {}

###################################################################################################

##	function trypsinDigest
##	input: protein sequence (string, case-insensitive)
##	output: tryptic peptides from input protein sequence (list of peptides, upper-case)
############# need to account for presence of proline #############
############# can probably filter out peptides below a certain length (maybe... 6 amino acids?) #############
def trypsinDigest(protSeq):
#	print("Returning peptides for given protein sequence ")		##remove
#	print(protSeq.upper())	##remove
	
	trypticPeptides = sampleProtein.replace('K','K_').replace('R','R_').split('_')

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
	protSeq = fastaRecord[1]

	##	parsing infoLine...

	print "InfoLine in Sample Fasta Record:"
	print infoLine

	print "Protein Sequence in Sample Fasta Record:"
	print protSeq
	##	...for NCBI Protein ID

	##	...for protein name

	##	...for organism ID
	splitInfoLine1 = infoLine.split('[')
	if ']' in splitInfoLine1[len(splitInfoLine1)-1]:	
		organismID = splitInfoLine1[len(splitInfoLine1)-1].strip(']')	## working for now, but need to account for entries without brackets in infoLine []
	else:
		organismID = None

	print "Organism ID: "
	print organismID		## remove
	## construct dictionary

	return parsedFastaRecord

###################################################################################################

parseFastaRecord(sampleFastaRecord)

###################################################################################################

##	function registerPeptides
##	input: FASTA record (modified; Python list)
##	input list format: [{'NCBIID':VALUE, 'proteinName':VALUE, 'organismID':VALUE, 'protSeq':VALUE}, [PEPTIDE1, PEPTIDE2, PEPTIDE3, ...]]
############# need to see what data should be extracted from the first > line of the FASTA entry... #############
############# eventually need to have more information about each key (more than just parent protein) -- start site in protein sequence, bacterial species/organism ID, etc. #############
def registerPeptides(fastaPeptides):
	for peptide in fastaPeptides[1]:
		if peptide not in allPeptides:

				allPeptides[peptide] = fastaPeptides[0]
		else:
			print "1 redundant"

	return allPeptides

###################################################################################################

print "allPeptides dictionary: "
print allPeptides
print "length of dictionary: "
print len(allPeptides)

registerPeptides(sampleFastaPeptides)

#import sys

#if len(sys.argv) != 2:
#	print 'Need one argument -- a filename'
#	print 'Proper usage: >>fastaparse.py somefile.fasta'
#	sys.exit()

#fileName = sys.argv[1]
#fileDict = {}
#inputFile = open(fileName)

#for fileLine in inputFile:
#	if '>' in inputFile.readline():
		

#inputFile.close()