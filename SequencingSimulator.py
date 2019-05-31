from scipy import stats
import random
import os
import numpy
import numbers

nucleotids = ['A', 'C', 'G', 'T']
complNucleotids = {
  "A": "T",
  "T": "A",
  "G": "C", 
  "C": "G"
}

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
    return "".join(read)

# Genetares 2 FASTQ files which contain paired-end reads.
def generateReads(refGenomeDict, quality, coverage, readSize, insertSize, fileName):
    samFile = open(fileName + ".sam", "w")
    for refGenomeName, refGenome in refGenomeDict.items():
        samFile.write("@SQ SN:{} LN:{}\n".format(refGenomeName, len(refGenome)))
    read1File = open(fileName + "_read1.fastq","w")
    read2File = open(fileName + "_read2.fastq", "w")
    for refGenomeName, refGenome in refGenomeDict.items():
        genomeSize = len(refGenome)
        fragmentNumber = numOfFragments(coverage, genomeSize, readSize)
        for i in range(fragmentNumber):
            # Insert size follows a normal distribution so we are simulating that.
            recalInsertSize = int(numpy.random.normal(insertSize, insertSize * 0.05))
            # Randomly choose position for the next fragment.
            fragmentPosition = genomeSize
            while (fragmentPosition + recalInsertSize >= genomeSize):
                fragmentPosition = int(random.random()*(genomeSize - recalInsertSize + 1))
        
            readId = ("%s_%d_%d_0:0_0:0_%d" %(refGenomeName,fragmentPosition, (fragmentPosition + recalInsertSize), i))
            
            # Generate first paried-end read.
            refPos = fragmentPosition
            readOne = generateSingleRead(readSize, refGenome, refPos, 1, lambda x: x)
            qualityOfRead = generateQuality(quality, readSize)

            read1File.write("@{}/1\n{}\n+\n{}\n".format(readId, readOne, qualityOfRead))
            samFile.write("{} {} {} {}\n".format(readId, refPos, readOne, qualityOfRead))
            
            # Generate second paried-end read.
            refPos = fragmentPosition + recalInsertSize
            readTwo = generateSingleRead(readSize, refGenome, refPos, -1, lambda x: complNucleotids[x])
            qualityOfRead = generateQuality(quality, readSize)

            read2File.write("@{}/2\n{}\n+\n{}\n".format(readId, readTwo, qualityOfRead))         
            samFile.write("{} {} {} {}\n".format(readId, refPos, readTwo, qualityOfRead))


    read1File.close()
    read2File.close()
    samFile.close()

# Generates mutations by inserting them in referenceGenome.
def insertMutations(refGenomeDict, errorSNV, errorInDel):
    for refGenome in refGenomeDict.values():
        for i in range(len(refGenome)):
            errorProb = random.random()
            if (errorProb <= (errorSNV + errorInDel)):
                if(errorSNV > 0 and errorProb <= errorSNV): # Simulate SNV.
                    # Randomly chooses new nucleotid and writes (index in nucleotied list)|(SNV mask) to refGenome. 
                    # Index value is bitwise OR-ed with mask for SNV so we can distingush it from insertion.
                    refGenome[i] = ((nucleotids.index(refGenome[i]) + random.randint(1,3))%4) | SNV 
                elif(errorInDel > 0): #Simulate INDEL
                    if (random.random() < 0.5): #Simulate insertion
                        # Randomly choose position from nucleotid list and add it in refGenome.
                        refGenome.insert(i, random.randrange(4))
                    else: # Simulate deletion.
                        # Replace nucleotid with mask for deletion.
                        refGenome[i] = DEL

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
            refGenomeName = line[1:].rstrip('\n')
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
    isValid = validateParameters(quality, coverage, readSize, insertSize, errorSNV, errorInDel)
    if(isValid):
        try:
            fileNameBase = os.path.basename(refGenomeFile)
            fileName = os.path.splitext(fileNameBase)[0]
            refGenomeDict = readGenome(refGenomeFile)
            # Check if we need to simulate errors.
            if(errorSNV + errorInDel > 0):
                insertMutations(refGenomeDict, errorSNV, errorInDel)
            generateReads(refGenomeDict, quality, coverage, readSize, insertSize, fileName)
        except FileNotFoundError:
            print("File not found. Please check the path and the name and try again.")


simulatePairedEndSequencing("proba.fasta", 70, 3, 7, 12)