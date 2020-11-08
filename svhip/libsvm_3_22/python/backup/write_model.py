import sys 
import multiprocessing
from decimal import Decimal
from svmutil import *
import math
from svm import *
from svm import __all__ as svm_all
from read_input import *
import scaling_config
import itertools


################################# Default parameters ###################
C_base = 2
Gamma_Base = 2

def_c_min = 0
def_c_high = 10
def_Gamma_min = 0
def_Gamma_high = 10
def_c_num = 20
def_g_num = 20
def_Gamma_Base = 2 
def_C_base = 2 
kernel_type = 2 #radial basis function in libsvm param
svm_class = 0   #c_svc in libsvm param
def_n_fold = 10 #iterations for cross-validation
def_nu = 0.5
def_outf = 'decision.model'

def_opt = [None, def_c_min, def_c_high, def_Gamma_min, def_Gamma_high, def_c_num, def_g_num, def_n_fold, def_nu, def_outf]

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
NOT FOR USE RIGHT NOW.

Alternative to main function, not fully implemented. Main use is giving a second
entrypoint for data tupels that were not parsed using the testset_creation
pipeline. 

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
        s = (u-l)/(num-1)
        vals = [float(l+i*s) for i in range(num)]
        return vals

    cs = values(c_low, c_up, num_c)
    gs = values(g_low, g_up, num_g)
    
    return sorted(itertools.product(cs,gs))
   
def parameters_grid_search( categories, values, parameters):
	
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
	process_count = multiprocessing.cpu_count()
	
	with multiprocessing.Pool(process_count) as process_pool:

		res_list = []
		argx = []
		for point in search_grid:
			c_val, g_val = point
			argx.append([categories, values, n_folds, c_val, g_val, nu])
		res_list = process_pool.map(call_svm_crossvalidate, argx)
		
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
		raise TypeError("Please recheck input format.")
	return problem[0], problem[1]

def parse_problem_instance(category_list, input_data):
	if not isinstance(category_list, (list, tuple)) or not isinstance(input_data, (list, tuple)):
		print("Invalid input format.")
		return None
	elif len(category_list)!=len(input_data):
		print("Invalid input format.")
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
		raise e
	return problem_set
  
def call_svm_trainer(categories, values, param_dict, outfile):
	parameters = svm_parameter()
	try:
		parameters.svm_type = 0
		parameters.kernel_type = kernel_type
		#parameters.n = n_fold
		parameters.gamma = param_dict.get('log_gamma')
		parameters.cost = param_dict.get('log_c')
		parameters.nu = param_dict.get('nu')
		problem = svm_problem(categories, values)
		m = svm_train(problem, parameters)
		
		svm_save_model(outfile, m)
	except Exception as e:
		raise e
	
	return 1
	
def call_svm_crossvalidate(argx):
	categories, values, n_folds, c_val, g_val, nu = argx
	parameters = svm_parameter()
	prob = svm_problem(categories, values)
	
	try:
		parameters.svm_type = 0
		parameters.kernel_type = kernel_type
		parameters.n = n_folds
		parameters.gamma = gamma_val(g_val)
		parameters.cost = C_val(c_val)
		parameters.nu = nu
		parameters.cross_validation = True
		MSE = svm_train(prob, parameters)
	except Exception as e:
		raise e
	
	return MSE 
	
###############################Main#####################################

def main(argx=None):
	"""
	Options:
	options = [None, c_low, c_high, g_low, g_high, num_c, num_g, n_fold, nu, output_name]
	"""
	parameters = svm_parameter()
	options = [None]*10
	inputfile = None
	
	if argx==None:
		argx = sys.argv
	else:
		pass
	
	if '-h' in sys.argv:
		print_help()
		return None
	
	for a in range(1, len(sys.argv)-1):
		if sys.argv[a] == '-c_low':
			try:
				options[1] = float(sys.argv[a+1])
			except Exception as e:
				raise TypeError("-c_low has to be float.")
		elif sys.argv[a] == '-c_high':
			try:
				options[2] = float(sys.argv[a+1])
			except Exception as e:
				raise TypeError("-c_high has to be float.")
		elif sys.argv[a] == '-g_low':
			try:
				options[3] = float(sys.argv[a+1])
			except Exception as e:
				raise TypeError("-g_low has to be float.")		  
		elif sys.argv[a] == '-g_high':
			try:
				options[4] = float(sys.argv[a+1])
			except Exception as e:
				raise TypeError("-g_high has to be float.")
		elif sys.argv[a] == '-num_c':
			try:
				options[5] = float(sys.argv[a+1])
			except Exception as e:
				raise TypeError("-num_c has to be float.")
		elif sys.argv[a] == '-num_g':
			try:
				options[6] = float(sys.argv[a+1])
			except Exception as e:
				raise TypeError("-num_g has to be float.")
		elif sys.argv[a] == '-n_fold':
			try:
				options[7] = float(sys.argv[a+1])
			except Exception as e:
				raise TypeError("n_fold has to be float.")
		elif sys.argv[a] == '-nu':
			try:
				options[8] = float(sys.argv[a+1])
			except Exception as e:
				raise TypeError("nu has to be float.")		
		elif sys.argv[a] == '-outf':
			try:
				options[9] = str(sys.argv[a+1])
			except Exception as e:
				raise e
	'''
	Last command line argument has to be path of input file:
	'''
	if isinstance(sys.argv[len(sys.argv)-1], str):
		try:
			tes = open(str(sys.argv[len(sys.argv)-1]), 'r')
			tes.close()
			inputfile =sys.argv[len(sys.argv)-1]
		except (FileNotFoundError, IOError):
			raise TypeError('Invalid input filename - please recheck directory and file.')
	else:
		raise TypeError("Last argument has to be input file. Filename has to be of Type String.")
		
	'''
	If options for grid search are not specified, default is loaded:
	'''
	for n in range(0, len(options)):
		if options[n] == None:
			options[n] = def_opt[n]		
	'''
	Set scaling ranges for this dataset:
	'''		
	z_range = scaling_config.get_z_range(str(inputfile) + '.scal')
	SCI_range = scaling_config.get_SCI_range(str(inputfile) + '.scal')
	shannon_range = scaling_config.get_shannon_range(str(inputfile) + '.scal')
	ranges = [z_range, SCI_range, shannon_range]
	'''
	Attempt to create problem instance:
	'''
	try:
		categories, values = read_problem_data(inputfile)
	except Exception as e:
		raise e
		
	'''
	Values will be scaled:
	'''
	values = scale_prob_data(values, ranges)
	print(values)
	'''
	Create grid search and write model:
	'''
	smallest_error, param_dict = parameters_grid_search(categories, values, options[1:9])
	ret = call_svm_trainer(categories, values, param_dict, options[9])
	if ret == 1:
		print("Process succesfully finished.")
	else:
		raise ValueError("Something went wrong. Please check input format and paramter range.")

################################# Misc. ################################

def print_help():
	with open('write_model.help', 'r') as f:
		for line in f.readlines():
			print(str(line))
	return None

########################################################################
if __name__ == '__main__':
    main()
	
def write_m(argx = None):
	main(argx)
