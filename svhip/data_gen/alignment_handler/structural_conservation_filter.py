#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 14:04:54 2020

@author: christopher
"""

from window_handle import window_handle
import tree_edit_distance as dist
import os
from math import floor
from random import sample

def get_control_group(directory):
    ls = []
    for entry in os.listdir(directory):
        if 'random' in entry and 'aln_align' in entry:
            ls.append(window_handle(os.path.join(directory, entry), False))
    return ls

def get_native_group(directory):
    ls = []
    for entry in os.listdir(directory):
        if 'random' not in entry and 'aln_align' in entry:
            ls.append(window_handle(os.path.join(directory, entry), True))
    return ls

def get_distances(aln_list):
    distances = []
    for aln in aln_list:
        distances.append(aln.get_distance())
    return distances

def get_cutoff(distances, k_value):
    distances_arrange = sorted(distances)
    cutoff_index = floor(k_value * len(distances_arrange))
    
    if cutoff_index < len(distances_arrange) -1:
        return distances_arrange[cutoff_index]
    elif len(distances_arrange) == 0:
        return 0
    else:
        return distances_arrange[0]

def reset_filter(aln_list):
    for a in range(0, len(aln_list)):
        aln_list[a].unmark()
    pass

def k_value_filter(aln_native, aln_control, k = 0.2):
    distances = get_distances(aln_native)
    distance_control = get_distances(aln_control)
    
    cutoff = get_cutoff(distance_control, k)
    aln_native_filtered = []
    
    for a in range(0, len(aln_native)):
        if aln_native[a].get_distance() >= cutoff:
            aln_native[a].mark_filtered()
        else:
            aln_native_filtered.append(aln_native[a])

    len_control = 0
    for a in range(0, len(aln_native)):
        if aln_native[a].get_filtered() != 1:
            len_control += 1

    if len(aln_control) >= len_control:
        aln_control = sample(aln_control, len_control)

    return aln_native_filtered, aln_control 