''' The preparation pipeline for trainingssets now packaged up in one 
module'''

import os
import sys

import_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
sys.path.append(import_path)
this_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_directory, 'alignment_handler/'))
sys.path.append(os.path.join(this_directory, 'help/'))
sys.path.append(os.path.join(this_directory, 'parser/'))
print(sys.path)

from window_handle import window_handle
from alignment_handle import Alignment_handle
from structural_conservation_filter import k_value_filter

sys.path.append(os.path.abspath(os.path.join(import_path, 'defaults/')))
sys.path.append(os.path.abspath(os.path.join(import_path, 'scaling/python/')))
sys.path.append(os.path.abspath(os.path.join(import_path, 'currysoup/')))
sys.path.append(os.path.abspath(os.path.join(import_path, 'logger/')))
sys.path.append(os.path.abspath(os.path.join(import_path, 'statistics/')))

import default_params
#from currysoup import creation_soup
import logger
from plotter import *
from cleanup import clean_dir

#########################Default parameters#############################
"""
These are default parameters, for options not set manually on main() call.
Feel free to change as needed in defaults/default_params.py.

default_category = '1'
default_minident = 50
default_maxident = 98
default_simulateAlignment = False
default_sampleSize = 100
default_simulateSampleSize = 1
default_num_processes = multiprocessing.cpu_count()
default_mute = False
"""
defaults = default_params.get_testset_cr_defaults()


############################Core function###############################
"""
The hub through which all above functions are called. Core function is 
called by main() once input cmd line is parsed. For options please see
section 'main function call' below.
"""

def create_training_set(filename, param, outfilename, statlog=None, gen_neg = 0):
    
    ret_value = 1
    minid, maxid, n_proc, mute = param[1], param[2], param[6], param[7]
    
    if '-' in param[0] or gen_neg == 0:
        neg_set = True
    else:
        neg_set = False
        
    '''
    Handlers for input files are initialized, files preprocessed and realigned.
    '''    
    aln = Alignment_handle(filename, native = True)
    aln.remove_all_gaps()
    aln.ident_Filter(maxid, n_proc, mute)
    aln.realign_me()
    
    '''
    Spawn a set of up to 10 partial alignments per 2..15 sequences. These
    'alignment windows' later form the basis for feature vector calculation.
    '''
    aln.spawn_subalignments(minid, maxid, n_proc)
    
    '''
    Generate control sets for alignment windows. SISSIz is used to simulate
    alignments with identical dinucleotide content and gap patterns. 
    '''
    control_alignments = []
    
    for i in range(0, len(aln.windows)):
        cw = aln.windows[i].generate_control_alignment()
        if verify_file(cw) is False:
            continue
        control_alignments.append(cw)
        
    '''
    Validation of structural conservation. Tree edit distances to 
    approximate structural distances between sequences are used as
    a metric to estimate overall structural conservation in comparison
    to the control set.
    This way, a sufficient separation between positive and negative 
    training instances on a structural level is ensured. 
    
    ToDo: Add k-cutoff value.
    '''
    aln_native, aln_negative = k_value_filter(aln.windows, control_alignments)
    
    '''
    Calculate feature vectors and save them in a .dat file. 
    '''
    for i in range(0, len(aln_native)):
        aln_native[i].calculate_feature_vector()
    for i in range(0, len(aln_negative)):
        aln_negative[i].calculate_feature_vector()
    write_data(aln_native, outfilename)
    write_data(aln_negative, outfilename)

    return ret_value
    
###########################main function call###########################
"""
Main functional procedure - input is a Fasta or clustal alignment file obtained from Rfam, preferrably either
regular alignment if |n_seq| < 300, otherwise seed alignment. Other input 
parameters, refer to -h or --man page.
"""

def main(argx=None, inputfile=None, outfile=None):
    THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
    ret_val = 0
    
    '''
    Should we call for help? 
    '''
    if '-h' in argx or '--man' in argx:
        print_help(THIS_FOLDER)
        return None
    '''
    Handle folder structure:
    '''
    if '/' in inputfile and '.ini' in inputfile:
        in_path = inputfile[:inputfile.rindex('/')+1]
    else:
        in_path = inputfile
    '''
    Initialize base parameters:
    '''
    if argx == None:
        argx = sys.argv
    else:
        pass
    '''
    Init statlog:
    '''
    statlog = creation_statlog(inputfile.split('/')[-1])    
    '''
    Parse function arguments with currysoup module:
    '''
    parameters, calls = creation_soup(argx, inputfile, outfile, function_log)
    """
    On main function call, at first input parameters are parsed. Obviously.
    If parameters are not set manually, default values are loaded.
    """

    for i in range(0, len(parameters)):
        if parameters[i] == None:
            parameters[i] = defaults[i]
    if outfile == None:
        if '.fasta' in str(inputfile):
            outfile = inputfile.replace('.fasta', '') + '.dat'
    
    """
    Should negative sets be auto-generated for all instances? Initialize boolean here, we will need this variable for later.
    """
    auto_gen = False
    
    """
    If any specific function calls to steps in the pipeline are made, 
    the necessary functions are called now:
    """
    if parameters[7] == False:
        print('Initialization complete. Starting main program...')
    """
    Core function is called with set options:
    """
    if '-readall' in calls:
        '''
        Check if input is .ini file in proper directory - alternative way to input positive and negative sets:
        '''
        if '.ini' in str(inputfile):
        
            try:
                cc = 0
                with open(inputfile) as  inf:
                    for line in inf.readlines():
                        if '#' in line:
                            break
                        else:
                            cc += 1
                            lin_block = line.split(':')
                            if '-' in lin_block[0]:
                                parameters[0] = '-1'
                                if parameters[7] == False:
                                    print('Reading File: ' + str(in_path + str(lin_block[1]).replace('\n', '')))
                                ret_val = create_training_set(in_path + str(lin_block[1]).replace('\n', ''), parameters, outfile, gen_neg = 0)
                            else:
                                parameters[0] = '1'
                                if parameters[7] == False:
                                    print('Reading File: ' + str(in_path + str(lin_block[1]).replace('\n', '')))
                                ret_val = create_training_set(in_path + str(lin_block[1]).replace('\n', ''), parameters, outfile, gen_neg = 0)
                    if parameters[7] == False:
                        print("End of input file reached. Data sets observed: " + str(cc))
                    if ret_val == 1:
                        if parameters[7] == False:
                            print("Process complete. Output written to " + str(outfile))
                            print("Process complete. Scaling range written to " + str(outfile)+ ' .scal')
                        return True
                    else:
                        lg_message = "Something went wrong. Please check input format. Since this is an unspecified error, please inform developer."
                        function_log.write_warning(lg_message) 
                        raise ValueError(lg_message)
            except Exception as e:
                function_log.write_log(str(e))
                raise e

        elif os.path.isdir(in_path):
            '''
            If valid is neither a valid .ini file with a list of input alignments and their respective categories nor a single alignment file, then it has to be a folder -
            should negative sets be autogenerated for all? 
            '''
            worktree = []
            
            y = 'n'
            #y = input("Automatically generate negative training instances for all input files? y/n \n")

            if y == 'y':
                auto_gen = True
            '''
            Build a list of alignment files with categories and traverse it:
            '''
            for entry in os.listdir(in_path):
                #category, file_list = read_input_dir(os.path.join(in_path, entry)) 
                #worktree.append((category, file_list))
                if parameters[7] is False:
                    print("Found file: " + str(entry))
                    worktree.append(os.path.join(in_path, entry))
            print('End of directory reached. Data sets found: ' + str(len(worktree)))
            
            for f in worktree:
                if auto_gen is True:
                    create_training_set(f, parameters, outfile, gen_neg = 1)
                else:
                    #y = input("Auto-generate negative set for " + str(f) + "? y/n \n")
                    if y == 'y':
                        create_training_set(f, parameters, outfile, gen_neg = 1)
                    else:
                        create_training_set(f, parameters, outfile, gen_neg = 0)
            ret_val = 1

            '''
            for a in range(0, len(worktree)):
                if worktree[a][1] is not None:         
                    for n in range(0, len(worktree[a][1])):
                        if worktree[a][1][n] is not None:
                            n_files += 1
                            if parameters[7] is False:
                                print('Reading File: ' + str(worktree[a][1][n]))
                            parameters[0] = worktree[a][0]
                            ret_val = create_training_set(worktree[a][1][n], parameters, outfile, statlog)
                print('End of directory reached. Data sets found: ' + str(n_files))
            '''
        else:
            lg_message = "FATAL: Input is neither valid fasta file,  .ini file or directory."
            e = TypeError(lg_message)
            function_log.write_log(lg_message)
            raise e
    
    elif '-extract' in calls:
        tripel =extract_data(inputfile)
        print(parameters[0]+' 1:'+str(tripel[0])+' 2:'+str(tripel[1])+' 3:'+str(tripel[2]))
        return True
    
    else:
        if auto_gen is False:
            y = input("Auto-generate negative set? y/n \n")
        if y == 'y' or auto_gen is True:
            ret_val = create_training_set(inputfile, parameters, outfile, gen_neg = 1)
        else:
            ret_val = create_training_set(inputfile, parameters, outfile, gen_neg = 0)
    '''
    The return value of create_training_set function specifies if an output
    file could be successfully written:
    '''
    if ret_val == 1:
        '''
        try:
            #statlog.draw_deletion_plot()
            #statlog.draw_distance_plot()
        except Exception:
            lg_message = "Could not write statistic log for graphic distance plot. Skipping..."
            function_log.write_warning(lg_message)
            print(lg_message)
        ''' 
        return True
    else:
        lg_message = "Something went wrong. Please check input format. Since this is an unspecified error, please inform developer."
        function_log.write_warning(lg_message) 
        raise ValueError(lg_message)
        
################################# Misc. ################################

def print_help(this_folder):
    with open(os.path.join(this_folder, 'help/generate_data.help'), 'r') as f:
        print(f.read())
    return None

def read_input_dir(dir_entry):
    seq_file_paths = []
    '''    
    If no category is given, ignore these folders:    
    '''    
    if os.path.isdir(dir_entry) and ('_neg_' in str(dir_entry) or '_pos_' in str(dir_entry)):
        for obj in os.listdir(dir_entry):
            if not os.path.isdir(os.path.join(dir_entry, obj)):
                seq_file_paths.append(os.path.abspath(os.path.join(dir_entry, obj)))
        '''
        Check if directory was empty:        
        '''    
        if len(seq_file_paths) == 0:
            return (None, None)
            '''
            Otherwise return filepaths of alignments and respective category:
            '''
        elif '_neg_' in str(dir_entry):
            return '-1',seq_file_paths
        elif '_pos_' in str(dir_entry):
            return '1',seq_file_paths
        '''
        Misc. error case/ should never be reached:
        '''
        raise IOError('Invalid folder - no category could be assigned.')
    else:
        return (None, None)
    
############################Outfile generation #########################
"""
Fills a .dat file with all generated feature vectors. This file can be used 
for classifier training.
"""

def write_data(windows, outfilename):
    
    if not outfilename.endswith('.dat'):
        outfilename = outfilename + '.dat'
    
    with open(outfilename, 'a') as dat:
        for i in range(0, len(windows)):
            values = windows[i].get_vector()
            if None in values or values == [0.0, 0.0, 0.0]:
                continue
            dat.write(str(windows[i].get_category()) +' 1:'+str(values[0])+' 2:'+str(values[1])+' 3:'+str(values[2])+'\n')
            
    return 1 
    
############################File verification###############################
"""
A helper function that checks if a file is empty. SISSIz mucks up sometimes.
"""

def verify_file(filename):

    if os.path.getsize(filename) == 0:
        lg_message = "WARNING: Could not create valid randomized alignment (empty file). Skipping, but rfam alignment should be rechecked."
        print(lg_message)
        #function_log.write_warning(lg_message)
        return False
    else:
        return True

########################################################################

#NOTE: Sorting is hierarchical from the top down i.e. first sequence is compared to all others, then second......
if __name__ == '__main__':
    main()
    
def testset_cr(argx = None, inputfile=None, outfile=None, flog=None):
    global function_log
    function_log = flog
    main(argx, inputfile, outfile)
    return None
