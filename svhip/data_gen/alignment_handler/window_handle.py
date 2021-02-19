
import data_gen.alignment_handler.tree_edit_distance as dist
import sys
from os import path
from os import pardir

import_path = path.abspath(path.join(path.dirname(path.abspath(__file__)), pardir))
sys.path.append(path.abspath(path.join(import_path, 'parser/')))
sys.path.append(path.abspath(path.join(import_path, 'RNAz_requirements/')))

from read_fasta import *
from data_gen.alignment_handler.RNAz_caller import extract_data
from subprocess import call
from shlex import split

''' 
This class handles the individual alignment windows generated from input alignments.
Originally it was the same class as the alignment_handle, but ultimately there
is no need for much of the function included in alignment_handle and
vice versa. 

The class itself is relatively light-weight and final, not much to change here.
'''

min_z, max_z = -8.15, 2.01
min_SCI, max_SCI = 0.0, 1.29
min_shannon, max_shannon = 0.0, 1.2878

class window_handle:
    
    identifier = ""
    path = ""
    native_group = False
    seqs = []
    tree_edit_distance = 0
    z_score = 0.0
    SCI = 0.0
    shannon = 0.0
    filtered = 0
    
    def __init__(self, path, native = True):
        self.identifier = path.split('/')[-1]
        self.seqs = parse_alignment_file(path)
        self.tree_edit_distance = dist.calculate_mean_distance(self.seqs)
        self.path = path
        self.filtered = 0
        if native is False:
            self.native_group = False
        else:
            self.native_group = True
            
    def calculate_feature_vector(self):
        values = extract_data(self.path) 
        self.z_score = self.scale_z(values[0])
        self.SCI = self.scale_SCI(values[1])
        self.shannon = self.scale_shannon(values[2])

    '''
      Caller Function for SISSIz: Generates a simulated alignment with comparable
      sequence composition. 
      Function parameters are the result of testing und should not be changed too much. 
      (Though increasing --flanks value can sometimes fix a crash on very large alignments
      - at the cost of run time.)
    '''

    def generate_control_alignment(self):
        outpath = self.path + ".random"
        with open(outpath, "w+") as outf:
            cmd = "SISSIz -n 1 -s --flanks 1750 " + self.path
            call(split(cmd), stdout= outf)
        return window_handle(outpath, native = False)


    '''
    Scaling is pretty straighforward - calculated values are normalized to a range of [-1  1] using the 
    minima/maxima extracted from the RNAz source code.
    
    Should the need arise, they may be changed at the top of the file.   
    '''
    
    def linear_scale(self, base, minimum, maximum):
        from_ = -1.0
        to_ = 1.0
        return from_+(to_ -from_)*(base -minimum)/(maximum - minimum)
    
    def scale_z(self, z_base):

        if z_base is not None:
            if z_base > max_z:
                z_base = max_z
            elif z_base < min_z:
                z_base = min_z
    
            return self.linear_scale(z_base, min_z, max_z)
    
        else:
            return 0.0

    def scale_SCI(self, SCI_base):
    
        if SCI_base is not None:
            if SCI_base > max_SCI:
                SCI_base = max_SCI
            elif SCI_base < min_SCI:
                SCI_base = min_SCI
    
            return self.linear_scale(SCI_base, min_SCI, max_SCI)
        else:
            return 0.0

    def scale_shannon(self, shannon_base):
    
        if shannon_base is not None:
            if shannon_base > max_shannon:
                shannon_base = max_shannon
            elif shannon_base < min_shannon:
                shannon_base = min_shannon
    
            return self.linear_scale(shannon_base, min_shannon, max_shannon)
        else:
            return 0.0

    def get_distance(self):
        return self.tree_edit_distance

    def get_vector(self):
        return [self.z_score, self.SCI, self.shannon]
    
    def get_category(self):
        if self.native_group == True:
            return str(1)
        else:
            return str(-1)
    
    def is_valid(self):
        if [self.z_score, self.SCI, self.shannon] != [0.0, 0.0, 0.0] and None not in  [self.z_score, self.SCI, self.shannon]:
            return True
        else:
            return False
    
    def get_filtered(self):
        return self.filtered
    
    def mark_filtered(self):
        self.filtered = 1

    def unmark(self):
        self.filtered = 0