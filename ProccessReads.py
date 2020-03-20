from Bio import SeqIO
import gzip
import os, sys
import subprocess
import regex
from gzip import open as gzopen
import numpy as np
import argparse
"""
Yogesh Goyal with Shweta Ramdas (started on March 23,2019) (most recent edit on June 13, 2019)
This program takes the two zip files of Read 1 and Read 2 from MiSeq sequencing and 
1. Looks for GFP forward primer seqnence (with an allowed specified error) in Read 2
2. Concatnates Read 1 and Read 2 correspinding to the ones that specify the condition #1. Creates a file "cellidsandbarcodes.txt".
3. Compares the CellID barcodes from 10X barcodes abd takes only those genes which match. Creates a file "matchingCellidsBarcodes.txt"


Modified by Michael Morley Added cmd line params 
    
"""

#Construct cmd line args 
# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-1", "--read1", required=True,
                help="1st Read Pair one with Cellid+UMI")
ap.add_argument("-2", "--read2", required=True,
                help="2nd Read Pair")
ap.add_argument("-s", "--staggerLength", required=True,type=int,default=1,
                help="Stagger Length")
ap.add_argument("-g", "--gfpPrimerLength",type=int,default=20,
                help="length of GFP primer")
ap.add_argument("-o", "--output", required=False,
                help="Name of the output directory")


args = ap.parse_args()

if args.output:
    ### Create output if needed
    try:
        os.mkdir(args.output)
    except OSError as error:
        print ("Creation of the directory %s failed" % args.output)
        print(error)
    else:
        print ("Successfully created the directory %s " % args.output)







staggerLength =1 #define the length of the stagger (e.g. For i7.3.1 corresponding to sample 1, it is 1)
gfpPrimerLength = 20 #including two extra bases AT after the GFP primer

read1 = []  #read 1, corresponding to CellID
read2 = []  #read 2, corresponding to GFP+ barcode
read1Qscore = [] #read 1 quality scores as array
read2Qscore = [] #read 2 quality scores as array

joinedRead1Read2 =[]
missingcellID = []
missingGfpBarcode = []
badQscore = []
screenedReads = []
shavedReads = []
badBarcode = []
testConditionID = []
uniqueScreenedReads = []
uniqueShavedReads = []
recordsR1 = SeqIO.parse(gzopen(args.read1, "rt"),format="fastq")
recordsR2 = SeqIO.parse(gzopen(args.read2, "rt"),format="fastq")

print("Matching reads...")
# Automatic check on whether the identifiers match for Read 1 and Read 2     
for record1, record2 in zip(recordsR1, recordsR2): 
        if record1.id == record2.id:
            read1.append(str(record1.seq))
            read2.append(str(record2.seq))
            read1Qscore.append(record1.letter_annotations["phred_quality"])
            read2Qscore.append(record2.letter_annotations["phred_quality"])

read2QscoreInt = np.array(read2Qscore) # convert it to array otherwise it gives error later
print("finished matching reads")            
joinedRead1Read2 = np.column_stack(((read1,read2))) # joining read 1 and read 2
print("joinRead1Read2Done")
print(len(joinedRead1Read2))
np.savetxt("joinCellIDsBarcodes.txt", joinedRead1Read2, fmt = '%s')

print("Starting filtering...")
for read, Qscore in zip(joinedRead1Read2, read2QscoreInt):
	read1Call = read[0]
	read2Call = read[1]
	cellID = read[0][0:16]
	UMI = read[0][17:28]
	barcode = read[1][21:91]  #needs to change depending on the sample
	shaved = np.array([cellID, UMI,barcode])
	read2QscoreCall = Qscore
	conditionGFP = 1
	conditionQscore = 1
	conditionBadBarcode = 1
     
	if len(regex.findall(r'(GGACGAGCTGTACAAGTAGG){e<=4}', read2Call[0:(args.staggerLength + args.gfpPrimerLength)])) == 0: # Do not consider reads without GFP tag
		missingGfpBarcode.append(read)
		conditionGFP = 0
    #Do not consider reads with phred Q score <14 for more than 5 times in GFP
	#if len(read2QscoreCall[np.where(read2QscoreCall[0:(staggerLength + gfpPrimerLength)] <= 14)]) > 5: 
	#	badQscore.append(read)
	#	conditionQscore = 0
    
	if (len(regex.findall("(AAAA)", read2Call)) > 0 or \
        len(regex.findall("(TTTT)", read2Call)) > 0 or \
        len(regex.findall("(GGGG)", read2Call)) > 0 or \
        len(regex.findall("(CCCC)", read2Call)) > 0 or \
        len(regex.findall("(NN)", read2Call)) > 0):
		badBarcode.append(read)
		conditionBadBarcode = 0
	if (conditionQscore == 1 and conditionGFP == 1 and conditionBadBarcode == 1):
		screenedReads.append(read)
		shavedReads.append(shaved)
        
uniqueScreenedReads = np.unique(screenedReads, axis = 0)
uniqueShavedReads = np.unique(shavedReads, axis = 0)
summaryFile = open(args.output + "/summaryFile.txt","w") 
summaryFile.write("total raw reads %d" % len(read1) + "\n")
#summaryFile.write("total missingcellID reads %d" % len(missingcellID)) 
summaryFile.write("total missingGfpBarcode reads %d" % len(missingGfpBarcode) + "\n") 
summaryFile.write("total badQscore reads %d" % len(badQscore) + "\n") 
summaryFile.write("total badBarcode reads %d" % len(badBarcode) + "\n") 
summaryFile.write("total screenedReads reads %d" % len(screenedReads) + "\n")
summaryFile.write("total screened Unique Reads reads %d" % len(uniqueScreenedReads) + "\n") 
summaryFile.write("total shavedReads reads %d" % len(shavedReads) + "\n")
summaryFile.write("total shaved Unique Reads reads %d" % len(uniqueShavedReads) + "\n") 
summaryFile.close()


#writing files to the Analysis folder
np.savetxt(args.output + "/joinedRead1Read2.txt", joinedRead1Read2, fmt='%s')
np.savetxt(args.output + "/missingcellID.txt", missingcellID, fmt='%s')        
#np.savetxt("testConditionID.txt", testConditionID)  
np.savetxt(args.output + "/badQscore.txt", badQscore, fmt='%s')  
np.savetxt(args.output + "/badBarcode.txt", badBarcode, fmt='%s')  
np.savetxt(args.output + "/screenedReads.txt", screenedReads, fmt='%s')
np.savetxt(args.output + "/uniqueScreenedReads.txt", uniqueScreenedReads, fmt='%s')  
np.savetxt(args.output + "/barcodesUniqueScreenedReads.txt", uniqueScreenedReads[:,1], fmt='%s')
np.savetxt(args.output + "/cellIDsUMIUniqueScreenedReads.txt", uniqueScreenedReads[:,0], fmt='%s')
np.savetxt(args.output + "/shavedReads.txt", shavedReads, fmt='%s')
np.savetxt(args.output + "/uniqueShavedReads.txt", uniqueShavedReads, fmt='%s')  
np.savetxt(args.output + "/barcodesUniqueShavedReads.txt", uniqueShavedReads[:,2], fmt='%s')
np.savetxt(args.output + "/UMIUniqueShavedReads.txt", uniqueShavedReads[:,1], fmt='%s')
np.savetxt(args.output + "/cellIDsUniqueShavedReads.txt", uniqueShavedReads[:,0], fmt='%s')
#unique UMIs

