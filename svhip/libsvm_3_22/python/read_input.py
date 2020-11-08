#reads and parses an RNAz input file for model generation
import sys
#sys.path.append('/libsvm-3.22/python/') ...get this to work, dammit.
from svmutil import *
import math
from svm import *
from svm import __all__ as svm_all


#################################Support subroutines####################
def isSequence(seq):
	if len(seq) > 1:
		nucleotides = set(['A','U','G','C','-','\n'])
		if set(seq) <= nucleotides:
			return True
		else:
			return False
	else:
		return False

def lookup_zScore(z):
	#THIS IS A PLACEHOLDER UNTIL THE ISSUE OF NON-NORMALDISTRIBUTED Z_SCORES IS SETTLED! A 
	#standard z-score table is not helping if actual distribution of scores is unknown and cannot
	#be reliably estimated.
	if abs(z) > 0:
		return z
	else:
		return float(z)
		
def math_log(n, base):
	if n != 0:
		return math.log(n, base)
	else:
		return 0

def calculate_ShannonEntropy(seqs):
	#Entropy is approximated by observing the total entropy of all individual columns using
	#the relative frequency of nucleotides per column. 
	N = 0.0
	N += len(seqs[0])
	svalues = []
	for i in range(0, len(str(seqs[0]))-1):
		nA, nU, nG, nC, gap, none = 0, 0, 0, 0, 0, 0
		nucs_dict = {
			"A": nA,
			"U": nU,
			"G": nG,
			"C": nC,
			"-": gap,
			"\n": none
		}
		#consensus sequence is ignored here
		for a in range(0, len(seqs)-1):
			key = seqs[a][i]
			literal_count = nucs_dict.get(key) + 1
			nucs_dict[key] = literal_count
		print([nucs_dict.get('A'), nucs_dict.get('U'), nucs_dict.get('G'), nucs_dict.get('C'),nucs_dict.get('-')])
		#l += sum([nucs_dict.get('A'), nucs_dict.get('U'), nucs_dict.get('G'), nucs_dict.get('C')])
		sA = float(nucs_dict.get('A')/N)*math_log(float(nucs_dict.get('A')/N),2)
		sU = float(nucs_dict.get('U')/N)*math_log(float(nucs_dict.get('U')/N),2)
		sG = float(nucs_dict.get('G')/N)*math_log(float(nucs_dict.get('G')/N),2)
		sC = float(nucs_dict.get('C')/N)*math_log(float(nucs_dict.get('C')/N),2)
		svalues.append(sum([sA, sG, sC, sU]))
		
	#print(svalues)
	return -float(1/N)* sum(svalues)

###############################read single RNAz output file#############

def read_rnaz_file(filename):
	#Takes an RNAz output file as input and extracts and normalizes z-score, SCI and 
	#Shannon-Entropy then returns them as a vector
	
	with open(filename) as f:
		i=0
		holdSeqs = []
		
		for line in f.readlines():
			if i==0:
				if 'Mean z-score:' in line:
					z = lookup_zScore(float(line.split(':')[1]))
					i+=1
			if i==1:
				if 'Structure conservation index:' in line:
					SCI = float(line.split(':')[1])
					i+=1
			else:
				if isSequence(str(line)):
					holdSeqs.append(line)
		SE = calculate_ShannonEntropy(holdSeqs)
	
	return [z, SCI, SE]
	
############################Read multiple rnaz files####################
#Expects a set of categories {-1, +1} and a corresponding set of filenames
#serving as training set. 

def read_multiple_rnaz(filenames):
	
	return_vector = []
	for i in range(0, len(filenames)):
		return_vector.append(read_rnaz_file(filenames[i]))

	return return_vector

#############################Compile SVM Problem file###################
#Just in case it is ever needed...

def compile_svm_problem(categories, values, filename):
	
	if len(categories)!=len(values):
		print("Mismatch in length of categories and values vector.")
		sys.exit(1)
	else:
		with open(filename, 'w') as f:
			for i in range(0, len(values)):
				f.write(str(categories[i]) + ' 1: ' +str(values[i][0])+ ' 2: ' +str(values[i][1])+ ' 3: ' +str(values[i][2]))
	
			

##########################Read regular SVM problem file#################

def call_svm_reader(filename):
	
	svm_problem_instance = svm_read_problem(filename)
	categories = svm_problem_instance[0]
	values = svm_problem_instance[1]
	return [categories, values]


###########TESTING######################################################
#a = read_rnaz_file('5SRNA.fa_block1.out')
#print(a)
########################################################################
