import subprocess
from SequencingSimulator import simulatePairedEndSequencing


def executeBwaMem(fileName):
    subprocess.run(["bwa", "index", "{}.fa".format(fileName)])
    subprocess.run(["bwa" ,"mem", "{}.fa".format(fileName), "{}_read1.fastq".format(fileName), "{}_read2.fastq".format(fileName)], stdout=open("{}_bwa.sam".format(fileName), "w"))

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
	values = [(0, 0), (0.01, 0), (0.05, 0), (0.1, 0), (0.2, 0), (0, 0.01), (0, 0.05), (0, 0.1), (0, 0.2), (0.01, 0.01), (0.05, 0.05), (0.15, 0.15)]
	resultsFile = open("results.txt", "w")
	for errorRate in values:
		fileName = "xxx"
		simulatePairedEndSequencing("{}.fa".format(fileName), 70, 4, 75, 100, errorRate[0], errorRate[1])
		executeBwaMem(fileName)
		result = compareSamFiles("{}.sam".format(fileName), "{}_bwa.sam".format(fileName))
		resultsFile.write("Eror rate for SNV: {}, error rate for INDEL: {}, BWA-MEM accuracy of alignment is: {}%.\n".format(errorRate[0], errorRate[1], result))

evaluate()
