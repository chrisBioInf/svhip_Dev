import sys 
import os
import_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))

if import_path not in sys.path:
    sys.path.append(import_path)
if 'libsvm_3_22/python' not in sys.path:
    sys.path.append(os.path.abspath(os.path.join(import_path, 'libsvm_3_22/python/')))
if 'logger/' not in sys.path:
    sys.path.append(os.path.abspath(os.path.join(import_path, 'logger/')))

import defaults.default_params as default_params
import multiprocessing
from decimal import Decimal
from svmutil import *
import math
from svm import *
from svm import __all__ as svm_all
from read_input import *
import itertools
from currysoup.currysoup import write_soup

################################# Default parameters ###################
'''
These initializes some parameters, for options not set manually on main() call.
Feel free to change as needed in defaults/default_params.py. 

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


##############################Parameter search grid#####################
"""
Creates a search grid of given ranges with intervals scaled according to 
desired number of steps (def: 20).
"""
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
    rel_pos = []
    d = (2**num_g) - 2

    for i in range(1, num_g +1):
        rel_pos.append((2**i )/d)

    def values(low, up, num, rel_pos):
        l = Decimal(str(low))
        u = Decimal(str(up))
        s = float(u-l)
        vals = [float((l)+Decimal(rel_pos[i]*s)) for i in range(int(num))]
        return vals

    cs = values(c_low, c_up, num_c, rel_pos)
    gs = values(g_low, g_up, num_g, rel_pos)
    

    return sorted(itertools.product(cs,gs)), cs, gs
   
def parameters_grid_search(categories, values, parameters, nproc, mute, epsilon, just_max = False, acc_name = "", recursive = False):
    
    c_low, c_high, g_low, g_high, num_c, num_g, n_folds, nu = parameters      
    
    search_grid, c_values, g_values = create_grid(c_low, c_high, g_low, g_high, num_c, num_g)
    #cn = (c_high - c_low) / (num_c - 1) / 2
    #gn = (g_high - g_low) / (num_g - 1) / 2
    process_count = nproc
    
    #Process the grid:
    with multiprocessing.Pool(process_count) as process_pool:

        res_list = []
        argx = []
        if mute is False:
            print("Begin " + str(n_folds) + "-fold cross validation...")
        
        for point in search_grid:
            c_val, g_val = point
            argx.append([categories, values, n_folds, c_val, g_val, nu, mute])
        res_list = process_pool.map(call_svm_crossvalidate, argx)
        print(res_list)
        
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
    
    '''
    as Felix correctly noted, we have a problem with errors being a hash table (due to likelyhood of duplicates). Extract maximum from
    res_list instead - implemented below

    As a side note, Siggi implemented a comparable algorithm at some point (but non-recursive). In his version,
    the hash-table-error issue remains and will probably lead to wrong results in practice. 
    Can fix this on demand, if deemed necessary.
    '''

    smallest_error = max(res_list)
    err_sorted = sorted(errors, reverse=True)
    
    print("Accuracy: " + str(smallest_error))
    best_c_val, best_g_val = errors.get(smallest_error) 
    
    '''
    If we do not want a recursive grid search, we can return the best parameter pair now. 
    Good for basic testing.
    '''
    if recursive == False:
        param_dict = {'log_c': int(best_c_val), 'log_gamma': int(best_g_val), 'nu': nu}
        return smallest_error, param_dict

    '''
    If this is a recursive check to estimate structure of the next grid, 
    only returning the absolute maximum is enough.
    '''
    if just_max is True:
        return max(res_list)
    
    if mute is False:
        print("Smallest error for: C = " + str(best_c_val) + ' gamma = ' + str(best_g_val))
    
    '''
    Helper functions to (1) define new flanks for the next recursive grid search and 
    (2) find the isolated maximum accuracy of the next grid. 
    '''
    def find_neighbors(point, values):

        low, high = point, point
        for i in range(0, len(values)):
            if values[i] == point:
                if (i-1) >= 0:
                    low = values[i-1]
                if (i+1) <= len(values)-1:
                    high = values[i+1]
        return low, high

    def find_deep_maximum(best_acc, res_list, argx, c_values, g_values, n_folds, nu, mute, nproc, categories, values, num_c, num_g, acc_name, epsilon, recursive):

        deep_accs = []
        cval = 0
        gval = 0
        n_screened = 0
        
        for n in range(0, len(res_list)):
            '''
            Multiple Maxima are very much possible and will be screened accordingly - but not more than 10 to somewhat limit run time.

            Proposal: The presence of already > 10 local maxima might already indicate a sufficient flattening of the accuracy fucntion topography - early exit point for runtime reduction?
            '''
            if res_list[n] == best_acc and n_screened <= 10:
                n_screened += 1
                cval, gval = argx[n][3], argx[n][4]
                
                c_high, c_low, g_high, g_low = new_flanks(cval, gval, c_values, g_values)
                parameters = [c_low, c_high, g_low, g_high, num_c, num_g, n_folds, nu]
                deep_accs.append(parameters_grid_search( categories, values, parameters, nproc, mute, epsilon, True, acc_name, recursive))
        
        '''
        Now, in case none of the checked sub-grids actually yield any increase 
        in accuracy, the process can be stopped then and there.
        '''
        if max(deep_accs) <= best_acc:
            return None
        
        deep_index = deep_accs.index(max(deep_accs))
        return deep_index                


    def new_flanks(best_c, best_g, c_values, g_values):
        #sorted_er = sorted(errors, reverse=True)
        lc, rc = find_neighbors(best_c, c_values)
        lg, rg = find_neighbors(best_g, g_values)
        #c_dist = abs(rg -lg)
        #g_dist = abs(rc -lc)

        #c_dist, g_dist are not really needed here, but where used in a test before
        return rc, lc, rg, lg

    '''
    Repeat grid search with closer values until we hit a set threshhold. Else return
    the found final values.
    
    Set epsilon to default if 0 or not specified (see /defaults/)
    '''
    if epsilon == 0.0:
        epsilon = default_params.get_epsilon()

    '''
    If errors are all equal, duplicate value to allow for dummy comparison in 
    the next step - otherwise too many redundant if-statements
    '''
    if len(err_sorted) == 1:
        err_sorted.append(err_sorted[0])
        
    '''
    Reiterate grid search with new flanks until maximal error difference of best values is
    below set epsilon:
    '''
    if abs(err_sorted[len(err_sorted) -1] - err_sorted[0]) > epsilon:
        deep_index = find_deep_maximum(max(res_list), res_list, argx, c_values, g_values, n_folds, nu, mute, nproc, categories, values, num_c, num_g, acc_name, epsilon, recursive)
        
        if deep_index is None:
            param_dict = {'log_c': int(best_c_val), 'log_gamma': int(best_g_val), 'nu': nu}
            return smallest_error, param_dict
        
        best_c_val, best_g_val = argx[deep_index][3], argx[deep_index][4]
        c_high, c_low, g_high, g_low = new_flanks(best_c_val, best_g_val, c_values, g_values)

        parameters = [c_low, c_high, g_low, g_high, num_c, num_g, n_folds, nu]
        return parameters_grid_search( categories, values, parameters, nproc, mute, epsilon, False, acc_name, recursive)
    
    else:
        param_dict = {'log_c': int(best_c_val), 'log_gamma': int(best_g_val), 'nu': nu}
        return smallest_error, param_dict

#########################Parse input data###############################
"""

First function is for processing data gathered with generate_data.py;
Second one is for manual compilation of categories and problem vectors
parsed with either read_input.py from RNAz files or concatenation of sets.
If possible, following the pipeline (first option) is advised.
"""
def read_problem_data(filename, function_log = None):
    try:
        problem = svm_read_problem(filename)
    except Exception as e:
        if function_log != None:
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
def create_svm_problem(categories, problem_values, function_log = None):
    try:
        problem_set = svm_problem(categories, problem_values)
    except Exception as e:
        if function_log is not None:
            function_log.write_log(str(e))
        raise e
    return problem_set
  
def call_svm_trainer(categories, values, param_dict, outfile):
    
    try:
        cmd = '-s {svm_type} -t {kernel_type} -c {c} -g {g} -b 1'.format(
            svm_type = 0,
            kernel_type = kernel_type,
            c = param_dict.get('log_c'),
            g = param_dict.get('log_gamma'),
            nu = param_dict.get('nu')
            )
        '''
        parameters.svm_type = 0
        parameters.kernel_type = kernel_type
        #parameters.n = n_fold
        parameters.gamma = param_dict.get('log_gamma')
        parameters.C = param_dict.get('log_c')
        parameters.nu = param_dict.get('nu')
        '''
        problem = svm_problem(categories, values)
        print('Svm trainer function of libsvm called...')
        m = svm_train(problem, cmd)
        if '.model' not in outfile:
            outfile = outfile + '.model'
        svm_save_model(outfile, m)
        print('File written: ' + str(outfile))
    except Exception as e:
        function_log.write_log(str(e))
        raise e
    
    return 1
    
def call_svm_crossvalidate(argx):
    categories, values, n_folds, c_val, g_val, nu, mute = argx
    
    cmd = '-s {svm_type} -t {kernel_type} -c {c} -g {g} -v {n}'.format(
            svm_type = 0,
            kernel_type = kernel_type,
            c = (c_val),
            g = (g_val),
            n = n_folds
            )
    
    if mute is False:
        print('Crossvalidation on datapoint: C = ' + str(c_val) + ', gamma = ' + str(g_val))
    '''
    NOTE: For completely unknown reasons the annotated code is defunct for now. 
    Implemented secondary solution that yields same
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
        ACC = svm_train(categories, values, cmd)
    
    except Exception as e:
        function_log.write_log(str(e))
        #print("Invalid gamma value... Skipping")
        return 0.0
    
    return ACC
    
###############################Main#####################################

def main(argx = None, inputfile=None, outfile=None):
    """
    options = [None, c_low, c_high, g_low, g_high, num_c, num_g, n_fold, nu, output_name]
    
    Initialize parameters, handle folder structure:
    """
    this_folder = os.path.dirname(os.path.abspath(__file__))
    parameters = svm_parameter()
    
    """
    Attempt to load function arguments if none are passed (i.e. direct function
    call was used. Deprecated, so this will not see much use...)
    """
    if argx == None:
        argx = sys.argv
    
    """
    Show help document for this specific module:
    """
    if '-h' in argx or '--man' in argx:
        print_help(this_folder)
        return None
    
    '''
    Parse function arguments with currysoup module:
    '''
    options, nproc, mute, epsilon, recursive_grid = write_soup(argx, inputfile, outfile, function_log)

    '''
    If options for grid search are not specified, default is loaded:
    '''
    for n in range(0, len(options)-1):
        if options[n] == None:
            options[n] = def_opt[n]    
    if nproc == None:
        nproc = multiprocessing.cpu_count()

    '''
    Attempt to parse training instance set:
    '''
    try:
        categories, values = read_problem_data(inputfile)
    except Exception as e:
        print("Failed to parse training data set. Is the path valid and are format specifications met?")
        #function_log.write_log(str(e))
        raise e
    '''
    Create grid search and write model:
    '''
    smallest_error, param_dict = parameters_grid_search(categories, values, options[1:9], 
                                                        nproc, mute, epsilon, False, 
                                                        inputfile.replace('.dat', '.acc'), recursive_grid)
    ret = call_svm_trainer(categories, values, param_dict, options[9])
    if ret == 1:
        print("Process succesfully finished.")
    else:
        lg_message = "Something went wrong. Please check input format and paramter range. Since this is an unspecified error, please inform developer."
        print(lg_message)
        #function_log.write_warning(lg_message) 
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
