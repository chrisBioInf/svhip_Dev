from Bio import SeqIO

'''
Quick scripts for parsing .fasta files and .aln files
'''

def read_fa(filename):
        try:
            seq_dict = {}
            records = SeqIO.parse(filename, "fasta")
            for record in records:
                seq_dict[str(record.id)] = str(record.seq)
            if len(seq_dict) >= 1:
                return seq_dict
            records = SeqIO.parse(filename, "clustal")
            for record in records:
                seq_dict[str(record.id)] = str(record.seq)
            return seq_dict
        except Exception:
            print("Invalid input file format detected. Neither fasta nor clustal alignment?")

def parse_alignment_file(filename, count = False):
    """
    Use with Bio.AlignIO
    """
    seq_list = []
    seq_n = int(list(reversed(filename.split('_')))[1])
    with open(filename) as align_f:
        alignment = SeqIO.parse(align_f, 'clustal')
        
        for record in alignment:
            if len(seq_list) < seq_n: 
                seq_list.append(str(record.seq))
            else:
                break
                
    if count == True:
        return len(seq_list)
    else:
        return seq_list

def isSequence(seq):
    if len(seq) > 1:
        nucleotides = set(['A','U','G','C', 'T', '-','\n'])
        if set(seq) <= nucleotides:
            return True
        else:
            return False
    else:
        return False
    
def count_seqs(filename):
    count  = 0
    with open(filename) as f:
        for l in f.readlines():
            if isSequence(l):
                count += 1
    return count