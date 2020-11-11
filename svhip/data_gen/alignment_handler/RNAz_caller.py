#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 09:58:29 2020

@author: christopher
"""

import multiprocessing
import shlex
import random
from math import log
from Bio import SeqIO
from re import sub 
from subprocess import call
from subprocess import check_output

##############################Draw Alignments############################
"""
For {2..15} sequences create 10 alignments with random average pairwise 
identity in interval d = [50..98].
Note: Max. deviation in seq. length should be 65%, even though that still 
seems pretty excessive / previous screening if necessary. --TODO

Note2 to self: Is it even realistically feasible to stay in the constraints
for sequence identity and have length deviation > 65%?
"""

def spawn_alignment(argx):
    alignment_file, avg_ident, n_seq, align_id, align_n, min_id = argx[0], argx[1], argx[2], argx[3], argx[4], argx[5]
    filename = alignment_file + '_align_' +  str(n_seq)+'_'+str(align_n)
    command_line = 'rnazWindow.pl --min-id=' + str(min_id) +' --opt-id='+str(avg_ident)+' -s 40 --no-reference --min-seqs='+str(n_seq)+' --max-gap=25 --max-seqs='+str(n_seq)+' '+str(alignment_file)
    arg_param = shlex.split(command_line)
    align_file = open(filename, 'w')
    call(arg_param, stdout = align_file)
    align_file.close()
    return [filename ,argx[3]]
    
def draw_alignments(alignment_file, min_ident, max_ident, num_processes, mute):
    align_dict = {}
    argx_list = []
    
    for n in range(2, 16):
        rand_idents = []
        for i in range(0, 10):
            avg_ident = random.randrange(min_ident, max_ident, 1)
            while avg_ident in rand_idents:
                avg_ident = random.randrange(min_ident, max_ident, 1)
            rand_idents.append(avg_ident)
            align_id = 'align_'+str(n)+'_'+str(i+1)
            seq_n = n
            align_n = i
            argx_list.append([alignment_file, avg_ident, seq_n, align_id, align_n +1, min_ident])
    
    with multiprocessing.Pool(num_processes) as pool:
        align_list = pool.map(spawn_alignment, argx_list)
    for entry in align_list:
        align_dict[entry[1]] = entry[0]
    if mute == False:
        print(str(len(align_dict))+ " Alignment blocks created.")        
    return align_dict

##############################calculate parameters######################
def math_log(n, base):
    if n != 0:
        return log(n, base)
    else:
        return 0

def calculate_ShannonEntropy(seqs):
    '''
    Entropy is approximated by observing the total entropy of all individual columns using
    the relative frequency of nucleotides per column. 
    '''
    N = 0.0
    N += len(seqs[0])
    svalues = []
    for i in range(0, len(str(seqs[0]))):
        nA, nU, nG, nC, gap, none = 0, 0, 0, 0, 0, 0
        nucs_dict = {
            "A": nA,
            "U": nU,
            "T": nU,
            "G": nG,
            "C": nC,
            "-": gap,
            "\n": none,
            " ": none
        }
        for a in range(0, len(seqs)):
            if i < len(seqs[a]):
                literal = nucs_dict.get(seqs[a][i])
                if literal is not None:
                    literal_count = literal + 1
                    nucs_dict[seqs[a][i]] = literal_count
                else:
                    pass
        #print([nucs_dict.get('A'), nucs_dict.get('U'), nucs_dict.get('G'), nucs_dict.get('C'),nucs_dict.get('-')])
        #l += sum([nucs_dict.get('A'), nucs_dict.get('U'), nucs_dict.get('G'), nucs_dict.get('C')])
        sA = float(nucs_dict.get('A')/N)*math_log(float(nucs_dict.get('A')/N),2)
        sU = float(nucs_dict.get('U')/N)*math_log(float(nucs_dict.get('U')/N),2)
        sG = float(nucs_dict.get('G')/N)*math_log(float(nucs_dict.get('G')/N),2)
        sC = float(nucs_dict.get('C')/N)*math_log(float(nucs_dict.get('C')/N),2)
        svalues.append(sum([sA, sG, sC, sU]))
        
    return -float(1/N)* sum(svalues)
    
def rnaz_call(filename):
    """
    This one is for now an ugly workaround. It effectively uses the RNAz
    routines to later extract MFE and z-Scores from the alignment without
    having to calculate them manually.
    Since RNAz uses a trained SVM to determine the z-scores we would 
    warp the results by using another (heuristic) algorithm.
    """
    command_line = 'RNAz ' + str(filename)
    byte_str = str(check_output(shlex.split(command_line)))
    
    rnaz_out = [None, None]
    rnaz_ls = byte_str.split(':')
    for i in range(0, len(rnaz_ls)):
        if rnaz_out[0] != None and rnaz_out[1] != None:
            break
        if "Mean z-score" in rnaz_ls[i]:
            rnaz_out[1] = truncate_bstr_element(rnaz_ls[i+1])
        if "Structure conservation index" in rnaz_ls[i]:
            rnaz_out[0] = truncate_bstr_element(rnaz_ls[i+1])

    return rnaz_out[0], rnaz_out[1]

def extract_data(filename):
    """
    Takes the sequences of exactly one alignment file and returns tripel
    of z-score, MFE, shannon-entropy.
    Note that these still have to be normalized later.
    """
    #print('Calculating values: ' + str(filename))
    seq_vector = parse_alignment_file(filename)
    if len(seq_vector) <= 1:
        print('Faulty alignment detected. Proceeding to next block.')
        return [None, None, None]
    
    SCI, z_score = rnaz_call(filename)
    shannon = calculate_ShannonEntropy(seq_vector)
    return [z_score, SCI, shannon]
    
def truncate_bstr_element(rnaz_entry):
    rnaz_out = sub('[a-zA-Z]','',rnaz_entry)
    return float(rnaz_out.replace("\\", ""))

##################Alignment helper #########################################

def parse_alignment_file(filename, count = False):
    """
    Use with Bio.AlignIO
    """
    seq_list = []
    #seq_n = int(list(reversed(filename.split('_')))[2])
    with open(filename) as align_f:
        alignment = SeqIO.parse(align_f, 'clustal')
        
        for record in alignment:
            seq_list.append(str(record.seq))                
    if count == True:
        return len(seq_list)
    else:
        return seq_list