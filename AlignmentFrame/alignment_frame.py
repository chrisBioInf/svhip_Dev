#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 07:42:06 2020

@author: christopher
"""

from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

#from RNAz_caller import *
#from window_handle import window_handle
from os import path
import numpy as np
import multiprocessing
from pandas import DataFrame
from pandas import Series
import matplotlib.pyplot as plt

class AlignmentFrame:
    
    aln_frame = None
    filepath = ""
    identifier = ""
    
    ###########################################################################
    """
    Core functions:
        Essential functions that handle the following:
        (a) first initialization of the frame from a sequence alignment file
        (b) update the frame, when it is expected that the file has changed
        (c) write changes from the frame into the managed file 
        (d) The parsing of the actual file I/O stream (with Biopython)
    """
    ###########################################################################
    
    def __init__(self, filepath, function_log =None):
        self.filepath = path.abspath(filepath)
        self.identifier = self.filepath.split('/')[-1]
        
        try:
            seq_dict = self.read_fa(self.filepath)
        except Exception as e:
            mssg = "Failed to initialize Alignment handle. Please recheck file. Path: " + str(path) + '\n'
            print(mssg)
            if function_log is not None:
                function_log.write_log(str(e))
                function_log.write_log(mssg)
            raise e
        
        seq_id, seq_list = list(seq_dict), []
        
        for key in seq_id:
            seq_list.append(seq_dict[key])
            
        self.aln_frame = DataFrame(data = {'Seq ID' : seq_id, 
                                       'Sequence' : seq_list})
        
        self.update_frame()
        
        
    def update_frame(self):
        self.add_sequence_lengths()
        self.add_CG_content()
        self.add_nucleotide_columns()
        self.add_gap_column()
        
        
    def update_file(self):
        with open(self.filepath, 'w+') as outf:
            for i in range(0, len(self.aln_frame['Sequence'])):
                seqid = self.aln_frame['Seq ID'].iloc[i]
                seq = self.aln_frame['Sequence'].iloc[i]
                outf.write('>' + seqid + '\n')
                outf.write(seq + '\n')
        #Most applications will require a manual call of realign_me()!        
        #self.realign_me()
        self.update_frame()
        
        
    def read_fa(self, filename):
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
    
    
    ###########################################################################
    """
    Whole Alignment manipulation:
        Interface for discarding sequences based on pairwise identity and 
        realigning the managed file with ClustalW2
    """
    ###########################################################################
        
    def identity_filter(self, maxident = 0.95, mute = False):
        
        discard_keys= []
        
        for i in range(0, len(self.aln_frame['Sequence']) -1):
            if mute == False:
                print('Scanning sequence index: ' + str(i))
            
            for a in range(i+1, len(self.aln_frame['Sequence'])):
                if calculateIdentity(self.aln_frame['Sequence'].iloc[i], self.aln_frame['Sequence'].iloc[a]) > maxident:
                    discard_keys.append(i)
                    break
                    
        self.aln_frame.drop(discard_keys, axis = 0)
        
        if mute == False:
            print("Total discarded sequences: " + str(len(discard_keys)))
            
        self.update_file()
        
        
    def realign_me(self):
        outpath = ""
        if self.filepath.endswith(".clw2"):
            outpath = self.filepath
        else:
            outpath = self.filepath + ".clw2"
        cline = ClustalwCommandline("clustalw2", infile = self.filepath, outfile = outpath)
        cline()
        self.filepath = outpath
        
        seq_dict = self.read_fa(self.filepath)
        seq_id, seq_list = list( seq_dict.keys() ), []
        for key in seq_id:
            seq_list.append(seq_dict[key])
        
        self.aln_frame['Sequence'] = seq_list
        self.aln_frame['Seq ID'] = seq_id
        self.update_frame()
        
    
    ###########################################################################
    """
    Sequence manipulation:
        
    """
    ###########################################################################

    def calculate_identity(self, str1, str2):
        
        def levenshtein(seq1, seq2):
            size_x = len(seq1) + 1
            size_y = len(seq2) + 1
            matrix = np.zeros ((size_x, size_y))
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
            return (matrix[size_x - 1, size_y - 1])
        
        ident = (1 - levenshtein(str1, str2)/max(len(str1), len(str2))) *100
        return ident
    
    def find(self, columnlabel, position):
        
        try:
            local_series = self.aln_frame[columnlabel]
            index = local_series[local_series == position].index[0]
            del(local_series)
            return index
        except Exception:
            e = ValueError("Unable to locate position " + str(position) 
            + " in column " + str(columnlabel) + "." )
            raise e
    
    def pairwise_identity(self, position = (None,None)):
        
        if type(position[0]) is int and type(position[1]) is int:
            seq1 = self.remove_gaps(self.aln_frame['Sequence'].iloc[position[0]])
            seq2 = self.remove_gaps(self.aln_frame['Sequence'].iloc[position[1]])
            
            return self.calculate_identity(seq1, seq2)
        
        if type(position[0]) is str and type(position[1]) is str:
            index1 = self.find('Seq ID', position[0])
            index2 = self.find('Seq ID', position[1])
            seq1 = self.remove_gaps(self.aln_frame['Sequence'].iloc[index1])
            seq2 = self.remove_gaps(self.aln_frame['Sequence'].iloc[index2])
            
            return self.calculate_identity(seq1, seq2)
        
        else:
            e = ValueError("Positional argument has to be given (either sequence index or ID)!")
            raise e 
            
    def drop_sequence(self, position = None):
        
        if type(position) == int:
            df = self.aln_frame.drop(position, axis = 0)
        
        elif type(position) == str:
            index = self.find('Seq ID', position)
            df = self.aln_frame.drop(index , axis = 0)
            
        else:
            e = ValueError("Positional argument has to be given (either sequence index or ID)!")
            raise e 
            
        self.aln_frame = df 
            
    def get_sequence(self, position = None, ungap = True):
        
        if type(position) is int:
            seq = self.aln_frame['Sequence'].iloc[position]
        
        elif type(position) is str:
            index = self.find('Seq ID', position)
            seq = self.aln_frame['Sequence'].iloc[index]
        
        else:
            e = ValueError("Positional argument has to be given (either sequence index or ID)!")
            raise e 
            
        if ungap is True:
            return self.remove_gaps(seq)
        else:
            return seq
    
    ###########################################################################
    """
    Counters, nucleotides and patterns:
        These functions evaluate and track nulceotide composition, seq length, 
        number of gaps, relative CG frequency and more.
        Important for statistics!
    """
    ###########################################################################
    
    def remove_gaps(self, sequence):
        return sequence.replace('-', '')
        
    #The following one discards all gaps in the actual Frame - not for use right now.
    def remove_all_gaps(self):
        local_series = []
        
        for i in range(0, len(self.aln_frame['Sequence'])):
            local_series.append(self.remove_gaps(self.aln_frame['Sequence'].iloc[i]))
        
        self.aln_frame['Sequence'] = np.asarray(local_series)
    
    def add_gap_column(self):
        gap_array = np.empty(len(self.aln_frame['Sequence']), dtype =int)
        
        for i in range(0, len(self.aln_frame['Sequence'])):
            gap_array[i] = self.aln_frame['Sequence'].iloc[i].count('-')
            
        self.aln_frame['gaps'] = gap_array
    
    def add_sequence_lengths(self):
        lengths = []
        
        for seq in self.aln_frame['Sequence']:
            lengths.append(len(self.remove_gaps(seq)))
            
        self.aln_frame['Length'] = lengths
        
    def add_A_column(self):
        array_A = np.empty(len(self.aln_frame['Sequence']), dtype =int)
        
        for i in range(0, len(self.aln_frame['Sequence'])):
            array_A[i] = self.aln_frame['Sequence'].iloc[i].count('A')
            
        self.aln_frame['A'] = array_A
        
    def add_C_column(self):
        array_C = np.empty(len(self.aln_frame['Sequence']), dtype =int)
        
        for i in range(0, len(self.aln_frame['Sequence'])):
            array_C[i] = self.aln_frame['Sequence'].iloc[i].count('C')
            
        self.aln_frame['C'] = array_C
        
    def add_G_column(self):
        array_G = np.empty(len(self.aln_frame['Sequence']), dtype =int)
        
        for i in range(0, len(self.aln_frame['Sequence'])):
            array_G[i] = self.aln_frame['Sequence'].iloc[i].count('G')
            
        self.aln_frame['G'] = array_G
        
    def add_U_column(self):
        array_U = np.empty(len(self.aln_frame['Sequence']), dtype =int)
        
        for i in range(0, len(self.aln_frame['Sequence'])):
            array_U[i] = self.aln_frame['Sequence'].iloc[i].count('U')
            
        self.aln_frame['U'] = array_U
        
    def add_nucleotide_columns(self):
        self.add_A_column()
        self.add_C_column()
        self.add_G_column()
        self.add_U_column()
        
    def get_dinucleotide_count(self, seq):
        if '-' not in seq:
            return float(len(seq) -1)
        else:
            return float(len(self.remove_gaps(seq)) -1)
        
    def add_CG_content(self):
        counter = 0
        cg_array = np.empty(len(self.aln_frame['Sequence']), dtype =float)
        
        for i in range(0, len(self.aln_frame['Sequence'])):
            seq_c = self.remove_gaps( self.aln_frame['Sequence'].iloc[i] )
            
            for b in range(0, len(seq_c) -1):
                if ( seq_c[b] + seq_c[b + 1] ) == 'CG':
                    counter += 1
            cg_array[i] = counter / self.get_dinucleotide_count(seq_c)
            counter = 0
        
        self.aln_frame['CG Content'] = cg_array
    
    def get_mean_length(self):
        return int(self.aln_frame['Length'].mean())
    
    def get_mean_cg(self):
        return self.aln_frame['CG Content'].mean()
    
    def get_mean_gaps(self):
        return self.aln_frame['gaps'].mean()
    
    def get_alignment_length(self):
        return len(self.aln_frame['Sequence'].iloc[0])
    
    def get_mean_A_abs(self):
        return int(self.aln_frame['A'].mean())
    
    def get_mean_C_abs(self):
        return int(self.aln_frame['C'].mean())
    
    def get_mean_G_abs(self):
        return int(self.aln_frame['G'].mean())
    
    def get_mean_U_abs(self):
        return int(self.aln_frame['U'].mean())
    
    def get_mean_A_freq(self):
        return float(self.aln_frame['A'].mean()) / float(self.get_mean_length())
    
    def get_mean_C_freq(self):
        return float(self.aln_frame['C'].mean()) / float(self.get_mean_length())
    
    def get_mean_G_freq(self):
        return float(self.aln_frame['G'].mean()) / float(self.get_mean_length())
    
    def get_mean_U_freq(self):
        return float(self.aln_frame['U'].mean()) / float(self.get_mean_length())
    
    ###########################################################################
    """
    Slicing:
        Returns different portions of the full frame as individual Series 
        objects and Slices. 
        - Useful for further processing elsewhere.
    """
    ###########################################################################
    
    def get_sequences(self):
        return self.aln_frame['Sequence']
    
    def get_seqid(self):
        return self.aln_frame['Seq ID']
    
    def get_lengths(self):
        return self.aln_frame['Length']
    
    def get_gaps(self):
        return self.aln_frame['gaps']
    
    def get_cg_series(self):
        return self.aln_frame['CG Content']
    
    def get_partial_frame(self, keywords):
        
        pf = DataFrame()
        
        for key in keywords:
            pf[key] = self.aln_frame[key]
        
        return pf
    
    def get_lighter_frame(self):
        return self.get_partial_frame(['Seq ID', 'Sequence', 'Length', 'gaps'])
    
    def subalignment(self, row_positions = (0, 0), column_positions = (0, 0), fp = None):
        
        if row_positions[1] == 0:
            row_positions[1] = len(self.aln_frame['Sequence'])
        
        if column_positions[1] == 0:
            column_positions[1] = len(self.aln_frame['Sequence'].iloc[0])
        
        if fp == None:
            fp= 'subset_' + str(row_positions[0]) + '_' + str(row_positions[1]) 
            + '_to_' + str(column_positions[0]) + '_' + str(column_positions[1])
            + '_' + self.filepath
        
        seq_list = list(self.aln_frame['Sequence'].iloc[row_positions[0] : row_positions[1]])
        seq_ids = list(self.aln_frame['Seq ID'].iloc[row_positions[0] : row_positions[1]])
        
        with open(fp, 'w+') as outf:
            for i in range(0, len(seq_list)):
                seqid = seq_ids[i]
                seq = self.remove_gaps(seq_list[i])
                seq = seq[column_positions[0], column_positions[1]]
                outf.write('>' + seqid + '\n')
                outf.write(seq + '\n')
        
        af = AlignmentFrame(fp)
        af.realign_me()
        return af
    
    def subset(self, row_1 = 0, row_2 = 0):
        
        if row_2 == 0:
            row_2 = len(self.aln_frame['Sequence'])
        
        df = self.aln_frame.copy()
        
        fp = 'subset_' + str(row_1) + '_' + str(row_2) + '_' + self.filepath
        df.filepath = fp
        
        if row_1 > 0:
            df.drop(np.arange(0, row_1, 1) , axis = 0)
        
        if row_2 < len(self.aln_frame['Sequence']):
            df.drop(np.arange(row_2, len(self.aln_frame['Sequence']), axis = 0))
        
        return df
    
    ###########################################################################
    """
    Plotting:
        Exactly what it says on the tin.
        Functions for plotting different axis of the frame - to be extended.
    """
    ###########################################################################
    
    def plot_length(self, save = False, name = None):
        #fig, ax = plt.subplots()
        
        y_high = max(self.aln_frame['Length']) + 5
        y_low = min(self.aln_frame['Length']) - 5
        
        plt.figure(figsize = (10, 8))
        ax = self.aln_frame.reset_index().plot(x = 'index', y = 'Length', kind = 'line',
                                  color = 'r', grid = True,
                                  ylim = (y_low, y_high),
                                  xticks = np.arange(0, len(self.aln_frame['Length']) + 5, 5))
        ax.set_xlabel('sequence index')
        ax.set_ylabel('sequence length [nt]')
        
        mean_l = self.get_mean_length()
        ax.plot([0, len(self.aln_frame['Length'])], [mean_l, mean_l],'k--', alpha = 0.7)
        
        if save == True and name != None:
            plt.savefig(name.replace('.pdf', '') + '.pdf', dpi = 200) 
        
        plt.show()
    
    
    def plot_CG(self, save = False, name = None):
        y_high = 1
        y_low = 0 
        
        plt.figure(figsize = (10, 8))
        ax = self.aln_frame.reset_index().plot(x = 'index', y = 'CG Content', kind = 'line',
                                  color = 'g', grid = True,
                                  ylim = (y_low, y_high),
                                  xticks = np.arange(0, len(self.aln_frame['CG Content']) + 5, 5))
        ax.set_xlabel('sequence index')
        ax.set_ylabel('frequency')
        
        if save == True and name != None:
            plt.savefig(name.replace('.pdf', '') + '.pdf', dpi = 200) 
        
        plt.show()
    
    
    def plot_gaps(self, save = False, name = None):
        y_high = max(self.aln_frame['gaps']) + 4
        y_low = min(self.aln_frame['gaps']) - 4
        
        plt.figure(figsize = (10, 8))
        ax = self.aln_frame.reset_index().plot(x = 'index', y = 'gaps', kind = 'line',
                                  color = 'b', grid = True,
                                  ylim = (y_low, y_high),
                                  xticks = np.arange(0, len(self.aln_frame['gaps']) + 5, 5))
        ax.set_xlabel('sequence index')
        ax.set_ylabel('gaps')
        
        mean_g = self.get_mean_gaps()
        ax.plot([0, len(self.aln_frame['Length'])], [mean_g, mean_g],'k--', alpha = 0.7)
        
        if save == True and name != None:
            plt.savefig(name.replace('.pdf', '') + '.pdf', dpi = 200) 
        
        plt.show()
    
    def plot_nucleotide_frequency(self, save = False, name = None):
        #fig, ax = plt.subplots()
        plt.figure(figsize = (10, 8))

        xs =  self.aln_frame.reset_index()['index']
        Ay = np.asarray(self.aln_frame['A'])
        Gy = np.asarray(self.aln_frame['G'])
        Cy = np.asarray(self.aln_frame['C'])
        Uy = np.asarray(self.aln_frame['U'])
        
        plt.plot(xs, Ay, 'b', alpha = 0.75, label = 'A')
        plt.plot(xs, Cy, 'g', alpha = 0.75, label = 'C')
        plt.plot(xs, Gy, 'r', alpha = 0.75, label = 'G')
        plt.plot(xs, Uy, 'm', alpha = 0.75, label = 'U')
        
        plt.xlabel('index')
        plt.ylabel('frequency')
        plt.xticks(np.arange(0, len(xs) + 5, 5))
        plt.legend(loc = 'upper right')
        
        if save == True and name != None:
            plt.savefig(name.replace('.pdf', '') + '.pdf', dpi = 200) 
        
        plt.show()
    
    def plot_nucleotide_hist(self, A = True, C = True, G = True, U = True, at_once = True, bins = 20):
        
        Ay = np.asarray(self.aln_frame['A'])
        Cy = np.asarray(self.aln_frame['C'])
        Gy = np.asarray(self.aln_frame['G'])
        Uy = np.asarray(self.aln_frame['U'])
        
        x_low = 0
        x_high = max((max(Ay), max(Cy), max(Gy), max(Uy)))
        
        if A is True:
            plt.figure(figsize = (8, 6))
            plt.hist(Ay, bins = np.arange(x_low, x_high + 10, 2), rwidth = 1, color = 'b', edgecolor = 'k')
            plt.xlabel('A [nt]')
            plt.ylabel('frequency')
            plt.xticks(np.arange(x_low, x_high + 10, 5))
            plt.show()
            plt.clf()
        
        if C is True:
            plt.figure(figsize = (8, 6))
            plt.hist(Cy, bins = np.arange(x_low, x_high + 10, 2), rwidth = 1, color = 'g', edgecolor = 'k')
            plt.xlabel('C [nt]')
            plt.ylabel('frequency')
            plt.xticks(np.arange(x_low, x_high + 10, 5))
            plt.show()
            plt.clf()
            
        if G is True:
            plt.figure(figsize = (8, 6))
            plt.hist(Gy, bins = np.arange(x_low, x_high + 10, 2), rwidth = 1, color = 'r', edgecolor = 'k')
            plt.xlabel('G [nt]')
            plt.ylabel('frequency')
            plt.xticks(np.arange(x_low, x_high + 10, 5))
            plt.show()
            plt.clf()
            
        if U is True:
            plt.figure(figsize = (8, 6))
            plt.hist(Uy, bins = np.arange(x_low, x_high + 10, 2), rwidth = 1, color = 'm', edgecolor = 'k')
            plt.xlabel('U [nt]')
            plt.ylabel('frequency')
            plt.xticks(np.arange(x_low, x_high + 10, 5))
            plt.show()
            plt.clf()
            
        if at_once is True:
            plt.figure(figsize = (8, 6))
            plt.hist(Ay, bins = np.arange(x_low, x_high + 10, 2), rwidth = 1, 
                     color = 'b', edgecolor = 'k', label = 'A', alpha = 0.5)
            plt.hist(Cy, bins = np.arange(x_low, x_high + 10, 2), rwidth = 1, 
                     color = 'g', edgecolor = 'k', label = 'C', alpha = 0.5)
            plt.hist(Gy, bins = np.arange(x_low, x_high + 10, 2), rwidth = 1, 
                     color = 'r', edgecolor = 'k', label = 'G', alpha = 0.5)
            plt.hist(Uy, bins = np.arange(x_low, x_high + 10, 2), rwidth = 1, 
                     color = 'm', edgecolor = 'k', label = 'U', alpha = 0.5)
            plt.legend(loc = 'upper left')
            plt.xlabel('nucleotides [n]')
            plt.ylabel('frequency')
            plt.xticks(np.arange(x_low, x_high + 10, 5))
            plt.show()
            plt.clf()
            
    def highlight_CG(self):
        self.highlight_subsequence('CG')
    
    def highlight_subsequence(self, subseq):
        
        seq_series = []
        id_series = []
        seq_array = []
        
        for a in range(0, len(self.aln_frame['Sequence'])):
            id_series.append(self.aln_frame['Seq ID'].iloc[a])
            seq_series.append(self.remove_gaps(self.aln_frame['Sequence'].iloc[a]))
            seq_array.append(seq_series[a].split( subseq ))
            
        seq_array = np.asarray(seq_array)
        
        for a in range(0, len(seq_array)):
            out_string = ""
            
            for c in range(0, len(seq_array[a]) -1):
                out_string = out_string + seq_array[a][c] + color.RED + subseq + color.END 
            
            out_string = out_string + seq_array[a][-1] 
            print(id_series[a] + ": \t" + out_string)
    
    ###########################################################################
    """
    Info functions:
        Helpful functions that print out diverse informations regarding the
        alignment frame and it's composition on the screen, such as 
        sequence and alignment length, nucleotide frequency and others, 
        including memory usage. 
    """
    ###########################################################################
    
    def info(self):
        self.aln_frame.info()
    
    def describe(self):
        self.aln_frame.describe()
        
    def print_me(self):
        print(self.aln_frame)
        
    def print_more(self):
        print(self.aln_frame)
        print("Mean sequence length: " + str(self.get_mean_length()))
        print("Mean CG frequency: " + str(self.get_mean_cg()))
        print("Length of alignment: " + str(self.get_alignment_length()))
        print("Total sequences managed: " + str(len(self.aln_frame['Sequence'])) + '\n')
        
        print("Average frequency of individual nucleotides (absolute/relative): ")
        print("A : " + str(self.get_mean_A_abs()) + ", " + str(round(self.get_mean_A_freq(), 4)) )
        print("G : " + str(self.get_mean_G_abs()) + ", " + str(round(self.get_mean_G_freq(), 4)) )
        print("C : " + str(self.get_mean_C_abs()) + ", " + str(round(self.get_mean_C_freq(), 4)) )
        print("U : " + str(self.get_mean_U_abs()) + ", " + str(round(self.get_mean_U_freq(), 4)) )
        
##############################################################################
"""
Helper functions for elimination of sequences with too high pairwise identity
and color formatting.
 Levenshtein distance is used as comparative metric for sequence identity. 
"""
##############################################################################

class color:
   PURPLE = '\033[1;35;48m'
   CYAN = '\033[1;36;48m'
   BOLD = '\033[1;37;48m'
   BLUE = '\033[1;34;48m'
   GREEN = '\033[1;32;48m'
   YELLOW = '\033[1;33;48m'
   RED = '\033[1;31;48m'
   BLACK = '\033[1;30;48m'
   UNDERLINE = '\033[4;37;48m'
   END = '\033[1;37;0m'


def levenshtein(seq1, seq2):
        size_x = len(seq1) + 1
        size_y = len(seq2) + 1
        matrix = np.zeros ((size_x, size_y))
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
    key_arr = np.array(seq_dict.keys())
    return np.array_split(key_arr, num_processes)
        
def spawn_search(argx):
    k, seq_dict, maxident, keylist_truncated = argx[0], argx[1], argx[2], argx[3]
    for key in keylist_truncated:
        if k != key and calculateIdentity(seq_dict.get(k), seq_dict.get(key)) > maxident:
            return k
        else:
            pass
    return None


#############################################################################
"""
TESTING:
"""
#############################################################################

af = AlignmentFrame('testing/test_1.fa.clw2')
af.print_more()


#af.highlight_CG()
#af.highlight_subsequence("GGGG")

#print(color.RED + "hello friends" + color.END)
'''
#####TEST SEQUENCE RETRIVAL:

seq = af.get_sequence(position = "X01003.1/1-119")
print(seq)

seq = af.get_sequence(position = 5)
print(seq)

af.drop_sequence(position = "X01003.1/1-119")
seq = af.get_sequence(position = 5)
print(seq)
'''

'''
labels = af.get_seqid()
print(labels)
'''
#af.ident_filter()
'''
### PLOTTING:

af.plot_length()
af.plot_CG()
af.plot_gaps()

af.plot_nucleotide_frequency()
af.plot_nucleotide_hist()
'''
#df = af.get_lighter_frame()

#print(af.aln_frame['Sequence'])
#print(af.aln_frame['Length'])
