#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 15:58:57 2020

A helper script that deletes all temporary alignment window files
--- currently not used.
"""
import os

def clean_dir(align_dict, mute):
    
    if mute is False:
        print("Beginning cleanup process for current run...")
    
    del_count = 0
    
    for key in align_dict:
        os.remove(align_dict.get(key))
        del_count += 1
        
    if mute is False:
        print("Current cleanup finished. " + str(del_count) +" temporary files deleted.")
    
    return None
