import sys 
import os
import_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
sys.path.append(import_path)
sys.path.append(os.path.abspath(os.path.join(import_path, 'libsvm_3_22/python/')))
sys.path.append(os.path.abspath(os.path.join(import_path, 'defaults/')))
sys.path.append(os.path.abspath(os.path.join(import_path, 'currysoup/')))
sys.path.append(os.path.abspath(os.path.join(import_path, 'logger/')))

import default_params
import multiprocessing
from decimal import Decimal
from svmutil import *
import math
from svm import *
from svm import __all__ as svm_all
from read_input import *
import scaling_config
import itertools
import subprocess
from currysoup import write_soup
import logger

################################# Default parameters ###################
'''
These are default parameters, for options not set manually on main() call.
Feel free to change as needed in defaults/default_params.py.

def_c_min = 0
def_c_high = 10
def_Gamma_min = 0
def_Gamma_high = 10
def_c_num = 20
def_g_num = 20
def_Gamma_Base = 2 
def_C_base = 2 
svm_class = 0   #c_svc in libsvm param
def_n_fold = 10 #iterations for cross-validation
def_nu = 0.5
def_outf = 'decision.model'
'''
C_base = 2
Gamma_Base = 2
kernel_type = 2 #radial basis function in libsvm param
def_opt = default_params.get_write_m_defaults()

def svm_type():
    return 'c_svc'

def gamma():
    #PLACEHOLDER - 
    return Gamma_Base
    
def c_exp(c):
    if c > 0:
        return math.log(c, C_base)
    else:
        return C_min

def gamma_exp(gamma):
    if gamma > 0:
        return log(gamma, Gamma_Base)
    else:
        return Gamma_min

def C_val(cexp):
    return C_base**cexp

def gamma_val(gexp):
    return Gamma_Base**gexp

"""
Function to set a few standard parameters, probably soon redundant.
"""
def set_parameters(svm_param):
    svm_param.gamma = gamma()
    svm_param.kernel_type = kernel_type
    svm_param.svm_type = svm_class
    return svm_param
    
################################## Scaling #############################
"""
SVM values taken into consideration are scaled to interval [-1, 1] - 
therefore the following functions have to be called before writing.
Modeled after scaling algorithm in RNAz.svm_helper.c i.e. training range
is assumed as absolute range.
No better solution in sight for now.
"""

def scale_z(z_base, z_range):

    return float(2*((z_base -z_range[0])/(z_range[1] -z_range[0])) -1)

def scale_SCI(SCI_base, SCI_range): 

    return float(2*((SCI_base -SCI_range[0])/(SCI_range[1] -SCI_range[0])) -1)
        
def scale_shannon(shannon_base, shannon_range):

    return float(2*((shannon_base -shannon_range[0])/(shannon_range[1] -shannon_range[0])) -1)
        
def scale_prob_data(dict_list, ranges):
    for dic in dict_list:
        dic[1]= scale_z(dic.get(1), ranges[0]) 
        dic[2]= scale_SCI(dic.get(2), ranges[1])
        dic[3]= scale_shannon(dic.get(3), ranges[2])
    return dict_list
        
def scale_vector(data_vector):
    hold_value = 0
    for i in range(0, len(data_vector)):
        hold_value = scale_z(data_vector[i][0])
        data_vector[i][0] = hold_value
        hold_value = scale_SCI(data_vector[i][1])
        data_vector[i][1] = hold_value
        hold_value = scale_shannon(data_vector[i][2])
        data_vector[i][2] = hold_value
    return data_vector

################################# Secondary function call ################
"""
NOT FOR USE RIGHT NOW. !!!!!!!!!!!!!!!!!!!

Alternative to main function, not fully implemented. Main use is giving a second
entrypoint for data tupels that were not parsed using the testset_creation
pipeline. 
With restructuring of package will probably be discarded soon.

#Each Data point has to be handed over as a 3-tuple of float (class, z-score, SCI, Shannon-Entropy). The
#Scores will have to be normalized. If parameters are to be set manually, they will have to be in dictionary form. The whole
#data vector has to be a list of these tuples.
"""
def write_model(categories, data_vector, name, parameters=None):

    if parameters== None:
        parameters_list = set_parameters(svm_parameter())
    elif isinstance(parameters, dict):
        parameters_list = svm_parameter()
        try:
            parameters_list.kernel_type = parameters.get('kernel_type')
            parameters_list.svm_type = parameters.get('svm_class')
            parameters_list.gamma = parameters.get('gamma')
            parameters_list.C = parameters.get('c_value')
        except KeyError:
            print('Paramters list given is incomplete.')
    else:
        print('Invalid parameter list given. Fallback to default.')
        parameters_list = set_parameters(svm_parameter())
    
    if not isinstance(data_vector, (list, tuple)):
        print("Invalid data input - cannot be parsed. Data vector has to be of type list.")
    else:
        for i in range(0, len(data_vector)):
            if not isinstance(data_vector[i], (list,tuple)):
                print("Invalid data input - cannot be parsed. Each Data point has to be handed over as a 3-tuple of float of format (z-score, SCI, Shannon-Entropy).")
                sys.exit(1)
            elif len(data_vector[i]) != 3:
                print("Invalid data input - cannot be parsed. Each Data point has to be handed over as a 3-tuple of float of format (z-score, SCI, Shannon-Entropy).")
                sys.exit(1)
            else:
                for a in range(0, 3):
                    if not isinstance(data_vector[i][a], float):
                        print('Invalid format in data tuple. Input values have to be float.')
                        sys.exit(1)
                    else:
                        continue
    
    #Scaling:
    data_scaled = scale_vector(data_vector)
    
    #############Call actual svm subroutines############################
    problem_vector = parse_problem_instance(categories, data_scaled)
    svm_problem_instance = create_svm_problem(problem_vector[0], problem_vector[1])
    #m = call_svm_trainer(svm_problem, parameters_list)
    m = call_svm_crossvalidate(svm_problem_instance, parameters_list, n_fold)
    
    svm_save_model(str(name), m)

##############################Parameter search grid#####################

def create_grid(c_low, c_up, g_low, g_up, num_c, num_g):
    """Create parameter search grid.

    Returns list of tuples. Each tuple has the form:
        (<log2(C) value>, <log2(gamma) value>)

    Arguments:
    - c_low -- Lower bound of log2(C).
    - c_up -- Upper bound of log2(C).
    - g_low -- Lower bound of log2(gamma).
    - g_up -- Upper bound of log2(gamma).
    - num_c -- Number of unique log2(C) values in the grid.
    - num_g -- Number of unique log2(gamma) values in the grid.
    """
    def values(low, up, num):
        l = Decimal(str(low))
        u = Decimal(str(up))
        s = ((u-l))/Decimal(num-1)
        vals = [float((l)+Decimal(i*s)) for i in range(int(num))]
        return vals

    cs = values(c_low, c_up, num_c)
    gs = values(g_low, g_up, num_g)
    
    return sorted(itertools.product(cs,gs))
   
def parameters_grid_search( categories, values, parameters, nproc, mute):
    
    c_low, c_high, g_low, g_high, num_c, num_g, n_folds, nu = parameters            
    '''
    Spawn several processes for optimal parameter estimation.
    Rightn now orients itself according to amount of cpu cores.
    Getting .apply_async from multiprocessing module to work would
    be desirable and will be further looked into.
    '''
    search_grid = create_grid(c_low, c_high, g_low, g_high, num_c, num_g)

    cn = (c_high - c_low) / (num_c - 1) / 2
    gn = (g_high - g_low) / (num_g - 1) / 2
    process_count = nproc
    
    with multiprocessing.Pool(process_count) as process_pool:

        res_list = []
        argx = []
        if mute is False:
            print("Begin " + str(n_folds) + "-fold cross validation...")
        
        for point in search_grid:
            c_val, g_val = point
            argx.append([categories, values, n_folds, c_val, g_val, nu, mute])
        res_list = process_pool.map(call_svm_crossvalidate, argx)
        
        if mute is False:
            print("Cross validation finished.")
        res_dict = {}
        for n in range(0, len(res_list)):
            point = argx[n][3], argx[n][4]
            res_dict[point] = res_list[n]
        errors = {}
        for point in search_grid:
            er = res_dict.get(point)
            errors[er] = point
    
    smallest_error = sorted(errors)[0]
    best_c_val, best_g_val = errors.get(smallest_error) 
    param_dict = {'log_c': best_c_val, 'log_gamma': best_g_val, 'nu': nu}
    
    if mute is False:
        print("Smallest error for: C = " + str(best_c_val) + ' gamma = ' + str(best_g_val))
    
    return smallest_error, param_dict

#########################Parse input data###############################
"""
First function is for processing data gathered with testset_creation.py;
Second one is for manual compilation of categories and problem vectors
parsed with either read_input.py from RNAz files or concatenation of sets.
If possible, following the pipeline (first option) is advised.
"""
def read_problem_data(filename):
    try:
        problem = svm_read_problem(filename)
    except Exception as e:
        function_log.write_log(str(e))
        raise e
    return problem[0], problem[1]

def parse_problem_instance(category_list, input_data):
    if not isinstance(category_list, (list, tuple)) or not isinstance(input_data, (list, tuple)):
        lg_message = "FATAL: Invalid format for input data. Expected are two lists/tuples with training categories and data vectors."
        function_log.write_warning(lg_message) 
        raise TypeError(lg_message)
        return None
    elif len(category_list)!=len(input_data):
        lg_message = "FATAL: Lengths of category list and data vector do not match!"
        function_log.write_warning(lg_message) 
        raise TypeError(lg_message)
        return None
    else:
        problem_vector = []
        problem_inner_vector = []
        problem_vector.append(category_list)
        for i in range(0, len(input_data)):
            new_dict = {
                        1: input_data[i][0],
                        2: input_data[i][1],
                        3: input_data[i][2]
                        }
            problem_inner_vector.append(new_dict)
        problem_vector.append(problem_inner_vector)
    
    return problem_vector

#########################Call libsvm subroutine#########################
def create_svm_problem(categories, problem_values):
    try:
        problem_set = svm_problem(categories, problem_values)
    except Exception as e:
        function_log.write_log(str(e))
        raise e
    return problem_set
  
def call_svm_trainer(categories, values, param_dict, outfile):
    parameters = svm_parameter()
    try:
        parameters.svm_type = 1
        parameters.kernel_type = kernel_type
        #parameters.n = n_fold
        parameters.gamma = param_dict.get('log_gamma')
        parameters.C = param_dict.get('log_c')
        parameters.nu = param_dict.get('nu')
        problem = svm_problem(categories, values)
        
        print('Svm trainer function of libsvm called...')
        m = svm_train(problem, parameters)
        
        svm_save_model(outfile, m)
        print('File written: ' + str(outfile))
    except Exception as e:
        function_log.write_log(str(e))
        raise e
    
    return 1
    
def call_svm_crossvalidate(argx):
    categories, values, n_folds, c_val, g_val, nu, mute = argx
    parameters = svm_parameter()
    prob = svm_problem(categories, values)
    
    cmd = '-s {svm_type} -t {kernel_type} -c {c} -g {g} -n {nu} -v {n}'.format(
            svm_type = 4,
            kernel_type = kernel_type,
            c = C_val(c_val),
            g = gamma_val(g_val),
            nu = nu,
            n = n_folds
            )
    
    if mute is False:
        print('Crossvalidation on datapoint: C = ' + str(c_val) + ', gamma = ' + str(g_val))
    '''
    NOTE: For completely unknown reasons the annotated code is defunct for now. 
    It worked before, I swear it. Implemented secondary solution that yields same
    results, so no real problem, but it remains odd.
    I will investigate when there is time.
    
    try:
        parameters.svm_type = 4
        parameters.kernel_type = kernel_type
        parameters.n = n_folds
        parameters.gamma = gamma_val(g_val)
        parameters.C = C_val(c_val)
        parameters.nu = nu
        parameters.cross_validation = True
        MSE = svm_train(prob, parameters)
    except Exception as e:
        print(str(e))
        function_log.write_log(str(e))
        raise e
    '''
    try:
        MSE = svm_train(categories, values, cmd)
    
    except Exception as e:
        function_log.write_log(str(e))
        raise e
    
    return MSE
    
###############################Main#####################################

def main(argx = None, inputfile=None, outfile=None):
    """
    options = [None, c_low, c_high, g_low, g_high, num_c, num_g, n_fold, nu, output_name]
    
    Initialize parameters, handle folder structure:
    """
    this_folder = os.path.dirname(os.path.abspath(__file__))
    parameters = svm_parameter()
    
    if argx == None:
        argx = sys.argv
    else:
        pass
    
    if '-h' in argx:
        print_help(this_folder)
        return None
    '''
    Parse function arguments with currysoup module:
    '''
    options, nproc, mute = write_soup(argx, inputfile, outfile, function_log)

    '''
    If options for grid search are not specified, default is loaded:
    '''
    for n in range(0, len(options)-1):
        if options[n] == None:
            options[n] = def_opt[n]    
    if nproc == None:
        nproc = multiprocessing.cpu_count()
    '''
    Set scaling ranges for this dataset:
    '''        
    if '.dat' in str(inputfile):
        z_range = scaling_config.get_z_range(str(inputfile).replace('.dat', '.scal'))
        SCI_range = scaling_config.get_SCI_range(str(inputfile).replace('.dat','.scal'))
        shannon_range = scaling_config.get_shannon_range(str(inputfile).replace('.dat','.scal'))
    else:
        z_range = scaling_config.get_z_range(str(inputfile) +'.scal')
        SCI_range = scaling_config.get_SCI_range(str(inputfile) + '.scal')
        shannon_range = scaling_config.get_shannon_range(str(inputfile) + '.scal')
    ranges = [z_range, SCI_range, shannon_range]
    '''
    Attempt to create problem instance:
    '''
    try:
        categories, values = read_problem_data(inputfile)
    except Exception as e:
        function_log.write_log(str(e))
        raise e
        
    '''
    Values will be scaled:
    '''
    values = scale_prob_data(values, ranges)
    '''
    Create grid search and write model:
    '''
    smallest_error, param_dict = parameters_grid_search(categories, values, options[1:9], nproc, mute)
    ret = call_svm_trainer(categories, values, param_dict, options[9])
    if ret == 1:
        print("Process succesfully finished.")
    else:
        lg_message = "Something went wrong. Please check input format and paramter range. Since this is an unspecified error, please inform developer."
        function_log.write_warning(lg_message) 
        raise ValueError(lg_message)

################################# Misc. ################################

def print_help(this_folder):
    with open(os.path.join(this_folder, 'write_model.help'), 'r') as f:
        print(f.read())
    return None

########################################################################
if __name__ == '__main__':
    main()
    
def write_m(argx = None, inputfile=None, outfile=None, flog=None):
    global function_log
    function_log = flog
    main(argx, inputfile, outfile)
    return None
