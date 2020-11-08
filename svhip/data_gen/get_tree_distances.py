import os
import sys 
import testset_creation as tc 
import statistics
import math
import random
import matplotlib.pyplot as plt

class window_handle:
    
    identifier = ""
    native_group = False
    seqs = []
    tree_edit_distance = 0
    z_score = 0.0
    SCI = 0.0
    shannon = 0.0
    filtered = 0
    
    def __init__(self, path, native):
        self.identifier = path.split('/')[-1]
        self.seqs = tc.parse_alignment_file(path)
        self.tree_edit_distance = tc.calc_dist_mean(self.seqs)
        '''
        values = tc.extract_data(path)
        self.z_score = self.scale_z(values[0])
        self.SCI = self.scale_SCI(values[1])
        self.shannon = self.scale_shannon(values[2])
        '''
        self.filtered = 0
        if native is False:
            self.native_group = False
        else:
            self.native_group = True
    
    def linear_scale(self, base, minimum, maximum):
        from_ = -1.0
        to_ = 1.0
        return from_+(to_ -from_)*(base -minimum)/(maximum - minimum)
    
    def scale_z(self, z_base):
        min_z, max_z = -8.15, 2.01
        if z_base is not None:
            if z_base > max_z:
                z_base = max_z
            elif z_base < min_z:
                z_base = min_z
    
            return self.linear_scale(z_base, min_z, max_z)
    
        else:
            return 0.0

    def scale_SCI(self, SCI_base): 
        min_SCI, max_SCI = 0.0, 1.29
    
        if SCI_base is not None:
            if SCI_base > max_SCI:
                SCI_base = max_SCI
            elif SCI_base < min_SCI:
                SCI_base = min_SCI
    
            return self.linear_scale(SCI_base, min_SCI, max_SCI)
        else:
            return 0.0

    def scale_shannon(self, shannon_base):
        min_shannon, max_shannon = 0.0, 1.2878
    
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
    
    def get_filtered(self):
        return self.filtered
    
    def mark_filtered(self):
        self.filtered = 1

    def unmark(self):
        self.filtered = 0

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
    cutoff_index = math.floor(k_value * len(distances_arrange))
    
    if cutoff_index < len(distances_arrange) -1:
        return distances_arrange[cutoff_index]
    elif len(distances_arrange) == 0:
        return 0
    else:
        return distances_arrange[0]

def append_vector(category, vector, outfile):
    if vector[0] == 0 and vector[1] == 0:
        return None

    with open(outfile, 'a') as f:
        if category == 1:
            f.write('-1 1:' +str(vector[0]) + ' 2:' + str(vector[1]) + ' 3:' + str(vector[2]) + '\n')
        else:
            f.write('1 1:' +str(vector[0]) + ' 2:' + str(vector[1]) + ' 3:' + str(vector[2]) + '\n')
    pass

def reset_filter(aln_list):
    for a in range(0, len(aln_list)):
        aln_list[a].unmark()
    pass

def k_value_filter(aln_native, aln_control, k, outfile):
    distance_control = get_distances(aln_control)
    cutoff = get_cutoff(distance_control, k)
    
    for a in range(0, len(aln_native)):
        if aln_native[a].get_distance() >= cutoff:
            aln_native[a].mark_filtered()
        else:
            append_vector(0, aln_native[a].get_vector(), outfile)

    len_control = 0
    for a in range(0, len(aln_native)):
        if aln_native[a].get_filtered() != 1:
            len_control += 1

    if len(aln_control) >= len_control:
        aln_control = random.sample(aln_control, len_control)

    for a in range(0, len(aln_control)):
        append_vector(1, aln_control[a].get_vector(), outfile)

    reset_filter(aln_native)

def main(path, k_arrange, outfile):
    #aln_native = get_native_group(path)
    #aln_native = get_native_group("non_filtered_control_group")
    #aln_control = get_control_group(path)
    aln_control = get_control_group("exp_1/set_1")
    
    for k in sorted(k_arrange):
        k_str = str(k)
        k_value_filter(aln_native, aln_control, k, outfile + '_k_' + k_str.replace('.', '') + '.dat' )

    print("Process finished.")

def draw_tree_edit_distances(path1, path2):
    aln1 = get_control_group(path1)
    aln2 = get_control_group(path2)
    tree1, tree2 = [], []
    for g in aln1:
        if g.get_distance() != 0:
            tree1.append(g.get_distance())
    for g in aln2:
        if g.get_distance() != 0:
            tree2.append(g.get_distance())
    
    fig, ax = plt.subplots()
    plt.hist(tree1, bins = 25, color = 'r', alpha = 0.5, edgecolor = 'k', label = "input control")
    plt.hist(tree2, bins = 25, color = 'b', alpha = 0.5, edgecolor = 'k', label = 'window control')
    plt.xlabel("mean tree \n edit distance")
    plt.ylabel("frequency")
    plt.subplots_adjust(bottom=0.15)
    plt.savefig("tree_dist.pdf")

draw_tree_edit_distances("exp_1/set_1","exp_2/set_1_control")
#main("", [0.1, 0.2, 0.3, 0.4, 0.5], "exp_1/control_group")
