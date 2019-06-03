import subprocess
from SequencingSimulator import simulatePairedEndSequencing


def executeBwaMem(fileName, errorRate):
    subprocess.run(["bwa", "index", "{}.fa".format(fileName)])
    subprocess.run(["bwa" ,"mem", "{}.fa".format(fileName), "{}_{}_read1.fastq".format(fileName, errorRate), "{}_{}_read2.fastq".format(fileName, errorRate)], stdout=open("{}_bwa.sam".format(fileName), "w"))

def compareSamFiles(ourSamFile, bwaSamFile):
	alignments = {}
	ourSam = open(ourSamFile, "r")
	correctlyAligned = 0;
	for line in ourSam.readlines():
		if line[0] == '@':
			continue;
		splitted = line.split('\t')
		mapId = splitted[0] + splitted[2]
		alignments[mapId] = splitted[1]
	ourSam.close()
	bwaSam = open(bwaSamFile, "r")
	for line in bwaSam.readlines():
		if line[0] == '@':
			continue;
		splitted = line.split('\t')
		mapId = splitted[0] + splitted[9]
		if alignments.get(mapId) == splitted[3]:
			correctlyAligned+=1
	return correctlyAligned/len(alignments)

def evaluate():
	values = [(0, 0), ()...] #add values you want to evaluate
	resultsFile = open("results.txt", "w")
	for errorRate in values:
		fileName = "" #add the name of refGenome file - without extension!
		simulatePairedEndSequencing("{}.fa".format(fileName), 70, 4, 150, 500, errorRate[0], errorRate[1])
		executeBwaMem(fileName, errorRate)
		result = compareSamFiles("{}_{}.sam".format(fileName, errorRate), "{}_bwa.sam".format(fileName))
		resultsFile.write("Eror rate for SNV: {}, error rate for INDEL: {}, BWA-MEM accuracy of alignment is: {}%.\n".format(errorRate[0], errorRate[1], result))

evaluate()
