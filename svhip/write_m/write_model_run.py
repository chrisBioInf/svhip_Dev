#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:29:07 2020

@author: christoktopus
"""

import write_model as wr
import sys

def main():
    opt = [None, '-g_low', 0.01, '-g_high', 16,'-c_low', 1, '-c_high', 4096, '-num_g',20,'-num_c', 20,'-n_fold' ,10]
    wr.write_m(opt, "sample_1.dat", "sample_1.model", None)
    
main()