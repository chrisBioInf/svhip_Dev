#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:29:07 2020

@author: christoktopus
"""

import write_model as wr
import sys

def main():
    opt = [None, '-c_low', 0.5, '-c_high',20,'-c_low', 0.5, '-c_high', 20, '-num_g',10,'-num_c', 10,'-n_fold' ,10]
    wr.write_m(opt, str(sys.argv[1]), str(sys.argv[1]) + '.model', None)
    
main()