from scipy import stats
import bisect
import random
import os
import numpy
import numbers
import time

nucleotids = ['A', 'C', 'G', 'T', 'N']
complNucleotids = {
  "A": "T",
  "T": "A",
  "G": "C", 
  "C": "G",
  "N": "N"
}

insertionPositions = {}

# Mask for SNV formated so the last 2 bits can be used for nucleotid index in nucleotids list
SNV = 0b100
DEL = 'D'

# Helper function for generating most optimal standard deviation for given mean.
# This is necessary because we are taking values from a speceific range of normal distribution.
def getStandardDeviation(averageQuality):
    if(averageQuality < 40 or averageQuality > 120):
        return 1
    elif(averageQuality < 50 or averageQuality > 105):
        return 2
    elif(averageQuality < 70 or averageQuality > 90):
        return 5
    else:
        return 10

# Generates qualities for nucleotids for one read using normal distribution and writes them in the given file 
def generateQuality(averageQuality, readSize):
    sigma = getStandardDeviation(averageQuality)
    a, b = (33 - averageQuality) / sigma, (126 - averageQuality) / sigma
    quality = stats.truncnorm.rvs(a, b, scale=sigma, loc=averageQuality, size=readSize)
    return "".join(map(lambda q: chr(int(round(q))), quality))

# Calculates necessary number of fragments based on the given coverage
def numOfFragments(coverage, genomeSize, readSize):
    return int(round((coverage * genomeSize) / (2 * readSize)))

# Finds the 1-based leftmost position where read matches refGenome.
def getLeftmostPosition(refGenome, startingPosition, endingPosition, refGenomeName):
    pos = startingPosition
    while((pos <= endingPosition) and isinstance(refGenome[pos], int) and  not(refGenome[pos] & SNV)):
        pos += 1
    return pos + 1 - bisect.bisect_right(insertionPositions.get(refGenomeName, []), pos);

# Generate one read. Direction defines should it be read from left to right or vice-versa.
def generateSingleRead(readSize, refGenome, refPos, direction, complFunction):
    j = 0;
    read = []
    while (len(read) < readSize):
        if(isinstance(refGenome[refPos], str)):
            if(refGenome[refPos] == DEL):
                refPos += direction
            else:
                read.append(complFunction(refGenome[refPos]))
                refPos += direction
        elif(refGenome[refPos] & SNV):
            read.append(complFunction(nucleotids[refGenome[refPos] & 3]))
            refPos += direction
        else:
            read.append(complFunction(nucleotids[refGenome[refPos] & 3]))
            refPos += direction
    return ("".join(read), refPos - direction)

# Genetares 2 FASTQ files which contain paired-end reads.
def generateReads(refGenomeDict, quality, coverage, readSize, insertSize, fileName, errorSNV, errorInDel):
    samFile = open("{}_({}, {}).sam".format(fileName, errorSNV, errorInDel), "w")
    for refGenomeName, refGenome in refGenomeDict.items():
        samFile.write("@SQ SN:{} LN:{}\n".format(refGenomeName, len(refGenome)))
        
    read1File = open("{}_({}, {})_read1.fastq".format(fileName, errorSNV, errorInDel),"w")
    read2File = open("{}_({}, {})_read2.fastq".format(fileName, errorSNV, errorInDel), "w")
    
    for refGenomeName, refGenome in refGenomeDict.items():
        print("[PREPROCESS]: Sequencing " + refGenomeName)
        genomeSize = len(refGenome)
        fragmentNumber = numOfFragments(coverage, genomeSize, readSize)
        for i in range(fragmentNumber):
            if i % 10000 == 0 and i > 0:
                print("[PREPROCESS]: created {} paired end reads".format(i))
            # Insert size follows a normal distribution so we are simulating that
            recalInsertSize = int(numpy.random.normal(insertSize, insertSize * 0.05))
            # Randomly choose position for the next fragment.
            fragmentPosition = genomeSize
            while (fragmentPosition + recalInsertSize >= genomeSize):
                fragmentPosition = int(random.random()*(genomeSize - recalInsertSize + 1))
            
            readId = ("{}_{}_{}_{}".format(refGenomeName, (fragmentPosition +1), (fragmentPosition + recalInsertSize + 1), i))
            
            # Generate first paried-end read.
            refPos = fragmentPosition
            readAndLastIndex = generateSingleRead(readSize, refGenome, refPos, 1, lambda x: x)
            qualityOfRead = generateQuality(quality, readSize)
            read1File.write("@{}/1\n{}\n+\n{}\n".format(readId, readAndLastIndex[0], qualityOfRead))

            leftmostPosition = getLeftmostPosition(refGenome, refPos, readAndLastIndex[1], refGenomeName)
            samFile.write("{}\t{}\t{}\t{}\n".format(readId, leftmostPosition, readAndLastIndex[0], qualityOfRead))
            
            # Generate second paried-end read.
            refPos = fragmentPosition + recalInsertSize
            readAndFirstIndex = generateSingleRead(readSize, refGenome, refPos, -1, lambda x: complNucleotids[x])
            qualityOfRead = generateQuality(quality, readSize)
            read2File.write("@{}/2\n{}\n+\n{}\n".format(readId, readAndFirstIndex[0], qualityOfRead)) 
            originalRead = []
            i = len(readAndFirstIndex[0]) - 1
            while(i > -1):
                originalRead.append(complNucleotids[readAndFirstIndex[0][i]])
                i-=1
            leftmostPosition = getLeftmostPosition(refGenome, readAndFirstIndex[1], refPos, refGenomeName)        
            samFile.write("{}\t{}\t{}\t{}\n".format(readId, leftmostPosition, "".join(originalRead), qualityOfRead))

    read1File.close()
    read2File.close()
    samFile.close()

# Generates mutations by inserting them in referenceGenome.
def insertMutations(refGenomeDict, errorSNV, errorInDel):
    for refGenomeName, refGenome in refGenomeDict.items():
        # Insert SNV.
        numOfSNV = round(len(refGenome)*errorSNV)
        while(numOfSNV > 0):
            errorPos = -1
            while(errorPos < 0 or isinstance(refGenome[errorPos], int)):
                errorPos = int(random.random()*len(refGenome))
            refGenome[errorPos] = ((nucleotids.index(refGenome[errorPos]) + random.randint(1,4))%5) | SNV
            numOfSNV -= 1
        
        # Insert INDEL.
        numOfInDel = round(len(refGenome)*errorInDel)
        print(numOfInDel)
        seqInsertions = []
        while(numOfInDel > 0):
            errorPos = -1
            while(errorPos < 0 or isinstance(refGenome[errorPos], int) or refGenome[errorPos] == 'D'):
                errorPos = int(random.random()*len(refGenome))
            if (random.random() < 0.5): #Simulate insertion
                 # Randomly choose position from nucleotid list and add it in refGenome.
                refGenome.insert(errorPos, random.randrange(5))
                seqInsertions.append(errorPos)
            else: # Simulate deletion.
                # Replace nucleotid with mask for deletion.
                refGenome[errorPos] = DEL
            numOfInDel -= 1
        seqInsertions.sort()
        insertionPositions[refGenomeName] = seqInsertions


# Loads genome from a file.
def readGenome(fileName):
    refGenome = []
    refFile = open(fileName)
    refGenomeName = ""
    refGenomeDict = {}
    for line in refFile.readlines():
        # Since one FASTA file can contain multiple sequences, we are loading them all and saving them in dict.
        if line[0] == '>':
            if refGenome:
                refGenomeDict[refGenomeName] = refGenome[:]
            refGenome = []
            refGenomeName = line[1:].rstrip('\n').split()[0]
        else:
            refGenome.extend(list(line.rstrip('\n')))
    if refGenome:
        refGenomeDict[refGenomeName] = refGenome[:]    
    refFile.close()
    return refGenomeDict
 
def checkPositiveIntValidity(value, name):
    if(not isinstance(value, int) or value < 0):
        print("{} must be whole number bigger then 0. Please choose the correct value and try again.".format(name))
        return 1
    return 0

def checkProbabilityValidity(value, name):
    if(not isinstance(value, numbers.Real) or value < 0 or value > 1):
        print("{} must represent probability between 0 and 1. Please choose the correct value and try again.".format(name))
        return 1
    return 0

def validateParameters(quality, coverage, readSize, insertSize, errorSNV, errorInDel):
    score = 0;
    if(not isinstance(quality, int) or quality < 26 or quality > 133):
        print("Quality must be whole number between 33 and 126. Please choose the correct value and try again")
        score = 1
    if(checkPositiveIntValidity(coverage, "Coverage") 
       | checkPositiveIntValidity(readSize, "ReadSize") 
       | checkPositiveIntValidity(insertSize, "InsertSize")):
        score = 1
    if(readSize > insertSize):
        print("ReadSize cannot be longer than InsertSize. Please choose the correct value and try again.")
        score = 1
    if(checkProbabilityValidity(errorSNV, "ErrorSNV") | checkProbabilityValidity(errorInDel, "ErrorINDEL")):
        score = 1
    if(errorSNV + errorInDel > 1):
        print("Sum of error rates for SNV  and INDEL cannot be larger than 1. Please choose the correct values and try again.")
        score = 1
    return score == 0
 
# Main simulation function.
def simulatePairedEndSequencing(refGenomeFile, quality, coverage, readSize, insertSize, errorSNV = 0, errorInDel = 0):
    print("Sequencing started for SNV: {}, InDel: {}".format(errorSNV, errorInDel))
    start = time.time()
    isValid = validateParameters(quality, coverage, readSize, insertSize, errorSNV, errorInDel)
    if(isValid):
        try:
            fileNameBase = os.path.basename(refGenomeFile)
            fileName = os.path.splitext(fileNameBase)[0]
            refGenomeDict = readGenome(refGenomeFile)
            # Check if we need to simulate errors.
            if(errorSNV + errorInDel > 0):
                insertMutations(refGenomeDict, errorSNV, errorInDel)
            generateReads(refGenomeDict, quality, coverage, readSize, insertSize, fileName, errorSNV, errorInDel)
        except FileNotFoundError:
            print("File not found. Please check the path and the name and try again.")
    print("Sequencing finished for SNV: {}, InDel: {}".format(errorSNV, errorInDel))
    end = time.time()
    print("Sequencing finished. Time elapsed: {}".format(end - start))