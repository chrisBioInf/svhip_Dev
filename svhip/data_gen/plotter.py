#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 16:33:30 2020

@author: christoktopus
"""

import matplotlib.pyplot as plt
import numpy as np
import statistics
import os

class creation_statlog:
    def __init__(self, filepath):
        self.family_deletion_dict = {}
        self.family_distance_dict = {}
        self.family_emptycount_dict = {}
        self.family_negative_distance_dict = {}
        self.filepath = "statistic_" + filepath
        try:
            os.mkdir(self.filepath)
        except Exception as e:
            print("Directory seems to exits. No biggie.")
    
    def deletion_log(self, filename, align_dict, rm_keys, mean_distances, mean_neg_distances, empty_files):
        with open("statistic/deleted_" + filename.split('/')[-1], 'w+') as fl:
            fl.write('Mean of structural distances over all blocks: ' +str(statistics.mean(mean_distances)) +'\n')
            if len(mean_neg_distances) < 1:
                mean_neg_distances = [0]

            fl.write('Mean of structural distances over all randomized blocks: ' +str(statistics.mean(mean_neg_distances)) +'\n')
            fl.write("Mean distances in order: " + str(mean_distances) +"\n")
            fl.write("Mean distances of randomized set: " +str(mean_neg_distances) +"\n")
            fl.write("Empty alignment blocks: " + str(empty_files) + "\n")
            fl.write('// \n')
            for key in rm_keys:
                fl.write(align_dict.get(key) + '\n')
            fl.write('// \n')
                
        self.family_deletion_dict[filename.split('/')[-1]] = len(rm_keys)
        self.family_distance_dict[filename.split('/')[-1]] = mean_distances
        self.family_negative_distance_dict[filename.split('/')[-1]] = mean_neg_distances
        self.family_emptycount_dict[filename.split('/')[-1]] = empty_files
        self.draw_hist(filename.split('/')[-1], mean_distances, mean_neg_distances)
        return None
        
    def draw_deletion_plot(self):
        fig, ax = plt.subplots()
        
        families = list(self.family_deletion_dict.keys())
        y, y2 = [], []

        for key in families:
            y.append(140 - (self.family_deletion_dict.get(key)+ self.family_emptycount_dict.get(key)))
            y2.append(140 -self.family_emptycount_dict.get(key))
        
        plt.ylim(0, 160)
        plt.yticks(ticks = np.arange(0, 160, 10))
        plt.title("Kept alignment blocks of Rfam families in dataset")
        plt.ylabel("Number of retained blocks")
        plt.xlabel("Name of Rfam familiy")
        plt.xticks(ticks = np.arange(0, len(families) +1, 1), labels = families)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        plt.plot(families, y, 'ro', alpha = 0.5, label="Kept alignment blocks")
        plt.plot(families, y2, 'b+', label="Total number of valid alignments")
        plt.grid()
        plt.legend()
        
        plt.savefig(self.filepath + "/deletion_plot.png", format = 'png',  dpi = 1600)
        
    def draw_distance_plot(self):
        fig, ax = plt.subplots()
        families = list(self.family_deletion_dict.keys())
        plt.title("Mean distances of all created alignment blocks")
        plt.ylabel("Mean distances \n[single edit operations]")
        plt.xlabel("Name of Rfam familiy")
        plt.xticks(ticks = np.arange(0, len(families) +1, 1), labels = families)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        
        for key in families:
            mean_distances = self.family_distance_dict.get(key)
            mean_neg_distances = self.family_negative_distance_dict.get(key)
            for mean in mean_distances:
                plt.plot(key, mean, 'r+', label='training data')
            for mean in mean_neg_distances:
                plt.plot(key, mean, 'b+', label='randomized alignment')
        
        plt.savefig(self.filepath +"/" +"distance_plot.png", format = 'png',  dpi = 1600)
        
    def draw_hist(self, family, mean_distances, mean_neg_distances):
        fig, ax = plt.subplots()
        
        plt.title("Mean structural distances comparison of \n positive sequence set and randomized alignment \n for family: " +str(family))
        plt.ylabel("frequency")
        plt.xlabel("mean distances \n [single edit operations]")
        
        bins = np.arange(0, 200 +1, 4)
        plt.xticks(ticks = np.arange(0, 200 +1, 10))
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        
        plt.hist(mean_distances, bins, color = 'b', label = 'training data', alpha = 0.5)
        plt.hist(mean_neg_distances, bins, color = 'r', label = 'randomized set', alpha = 0.5)
        plt.legend()
        plt.grid()
        
        plt.savefig(self.filepath +"/edit_distance_" +str(family) +".png", format = 'png',  dpi = 1600)
        return None
        
        
