from Bio import SeqIO
import sys

'''
Quick script for parsing .fasta files
'''

def read_fa(filename):

        try:
            seq_dict = {}
            records = SeqIO.parse(filename, "fasta")
            for record in records:
                seq_dict[str(record.id)] = str(record.seq)
            if len(seq_dict) >= 1:
                return seq_dict
            seq_dict = {}
            records = SeqIO.parse(filename, "clustal")
            for record in records:
                seq_dict[str(record.id)] = str(record.seq)
            return seq_dict
        except Exception:
            print("Invalid input file format detected. Neither fasta nor clustal alignment?")


