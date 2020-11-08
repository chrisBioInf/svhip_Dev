import sys
from Bio import SeqIO

def clustal_to_fa(filename, outf):
	records = SeqIO.parse(filename, "clustal")
	with open(outf, 'w') as outf:
		for record in records:
			outf.write('>' + str(record.id) + '\n')
			outf.write(str(record.seq) + '\n')
		'''
		count = SeqIO.write(records, filename + 'fasta', "fasta")
		print("Converted %i records" % cou
                '''
	return None

