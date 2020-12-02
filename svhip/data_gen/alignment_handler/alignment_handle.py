from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

from RNAz_caller import *
from window_handle import window_handle
import sys
from os import path
from os import pardir
import numpy as nump
import multiprocessing
from subprocess import call
from shlex import split

import_path = path.abspath(path.join(path.dirname(path.abspath(__file__)), pardir))
sys.path.append(path.abspath(path.join(import_path, 'parser/')))
from read_fasta import *

class Alignment_handle:
    
    identifier = ""
    path = ""
    native = True
    ids = []
    seq_dict = {}
    windows = []
    
    def __init__(self, path, native):
        self.seq_dict = read_fa(path)
        self.path = path
        #self.remove_all_gaps()
        self.identifier = path.split('/')[-1]
        self.native = native
    
    def remove_gaps(self, sequence):
        return sequence.replace('-', '')
        
    def remove_all_gaps(self):
        for key in self.seq_dict:
            tmp = self.remove_gaps(self.seq_dict[key])
            self.seq_dict[key] = tmp

    def ident_Filter(self, maxident, num_processes, mute):
        #key_matrix = split_keys(seq_dict, num_processes)
        discard_keys = []
        argx_list = []
        keylist = list(self.seq_dict.keys())
        
        for key in keylist:
            argx_list.append([key, self.seq_dict, maxident, keylist[keylist.index(key):]])
        
        with multiprocessing.Pool(num_processes) as process_pool:
            discard_keys = process_pool.map(spawn_search, argx_list)
        
        discard_n = 0
        for k in discard_keys:
            if k !=None:
                if mute == False:
                    print("Discarded sequence: " + self.seq_dict.get(k))
                discard_n += 1
                self.seq_dict.pop(k)
            else:
                pass
        if mute == False:
            print("Total discarded sequences: " + str(discard_n))
            
        self.update_file()
            
    def realign_me(self):
        outpath = ""
        if self.path.endswith(".clw2"):
            outpath = self.path
        else:
            outpath = self.path + ".clw2"
        
        cline = ClustalwCommandline("clustalw2", infile = self.path, outfile = outpath)
        cline()
        self.path = outpath
    
    def update_file(self):
        with open(self.path, 'w+') as outf:
            for key in self.seq_dict.keys():
                outf.write('>' + str(key) + '\n')
                outf.write(str(self.seq_dict[key]) + '\n')
        #Most applications will require a manual call of realign_me()!        
                
        #self.realign_me()
        
    def generate_control_alignment(self):
        outpath = self.path.replace(".clw2" , "") + ".random"
        with open(outpath, "w+") as outf:
            cmd = "SISSIz -n 1 -s --flanks 1750 " + self.path
            call(split(cmd), stdout= outf)
        return Alignment_handle(outpath, native = False)
    
    def spawn_subalignments(self, min_ident, max_ident, n_proc, window_size=120, step=40):
        align_dict = draw_alignments(self.path, min_ident, max_ident, n_proc, True, window_size, step)
        for key in align_dict.keys():
            self.windows.append(window_handle(align_dict[key], self.native))
            
    def return_windows(self):
        return self.windows

##############################################################################
"""
Helper functions for elimination of sequences with too high pairwise identity. Levenshten 
distance is used as comparative metric for sequence identity. 
"""
##############################################################################

def levenshtein(seq1, seq2):
        size_x = len(seq1) + 1
        size_y = len(seq2) + 1
        matrix = nump.zeros ((size_x, size_y))
        for x in range(0, size_x):
            matrix [x, 0] = x
        for y in range(0, size_y):
            matrix [0, y] = y
    
        for x in range(1, size_x):
            for y in range(1, size_y):
                if seq1[x-1] == seq2[y-1]:
                    matrix [x,y] = min(
                        matrix[x-1, y] + 1,
                        matrix[x-1, y-1],
                        matrix[x, y-1] + 1
                    )
                else:
                    matrix [x,y] = min(
                        matrix[x-1,y] + 1,
                        matrix[x-1,y-1] + 1,
                        matrix[x,y-1] + 1
                    )
        #print (matrix)
        return (matrix[size_x - 1, size_y - 1])

def calculateIdentity(str1, str2):
    ident = (1 - levenshtein(str1, str2)/max(len(str1), len(str2))) *100
    return ident
        
def split_keys(seq_dict, num_processes):
    key_arr = nump.array(seq_dict.keys())
    return nump.array_split(key_arr, num_processes)
        
def spawn_search(argx):
    k, seq_dict, maxident, keylist_truncated = argx[0], argx[1], argx[2], argx[3]
    for key in keylist_truncated:
        if k != key and calculateIdentity(seq_dict.get(k), seq_dict.get(key)) > maxident:
            return k
        else:
            pass
    return None
