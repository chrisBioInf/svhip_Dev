#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 11:56:25 2020

@author: christopher
"""

from generate_data import create_training_set

def main(inputfile, outputfile):
    
    parameters = ["1", 50, 98,0,0,0, 4, False ]
    ret_val = create_training_set(inputfile, parameters, outputfile)
    
main('testsets/rf00001.fasta', 'test.dat')