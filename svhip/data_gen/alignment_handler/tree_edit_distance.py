#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 13:56:32 2020

@author: christopher

A simple module that handles the calculation of a structural conservation approximation 
using a tree representation of the secondary structure. 
Depends on the ViennaRNA package. 
"""

import statistics
import RNA as vienna

def fold_rna(seq):
    fc = vienna.fold_compound(seq)
    mfe_struct, mfe = fc.mfe()
    seq_ex = vienna.expand_Full(mfe_struct)
    return [seq_ex, mfe]
    
def calculate_tree_edit_dist(seq1, seq2):
    seq_struc1, seq_struc2 = fold_rna(seq1)[0], fold_rna(seq2)[0]
    seq_tree1, seq_tree2 = vienna.make_tree(seq_struc1), vienna.make_tree(seq_struc2)
    return vienna.tree_edit_distance(seq_tree1, seq_tree2)

def calculate_mean_distance(seq_list):
    '''
    In case there are empty alignments created by RNAzWindow, which happens
    depending on set parameters for that run, 0 is returned. Zero-values are
    later removed from distance vector.
    Note: "Real" zero-alignments i.e. all sequences are equal can only be
    created if all input sequences in total are equal, so that should not be an
    issue.
    '''
    if len(seq_list) > 1: 
        distance_vector =[]
        for i in range(0, len(seq_list)-1):
            for a in range(i+1, len(seq_list)):    
                seq_struc1, seq_struc2 = fold_rna(seq_list[i])[0], fold_rna(seq_list[a])[0]
                seq_tree1, seq_tree2 = vienna.make_tree(seq_struc1), vienna.make_tree(seq_struc2)
                distance_vector.append(vienna.tree_edit_distance(seq_tree1, seq_tree2))
        return statistics.mean(distance_vector)
    else:
        return 0
