import multiprocessing

###########################Defaults: write_model.py#####################
def get_write_m_defaults():
	C_base = 2
	Gamma_Base = 2

	def_c_min = 0.03125
	def_c_high =65536
	def_Gamma_min = 0.00003051757
	def_Gamma_high = 32
	def_c_num = 20
	def_g_num = 20
	def_Gamma_Base = 2 
	def_C_base = 2 
	kernel_type = 2 #radial basis function in libsvm param
	svm_class = 0   #c_svc in libsvm param
	def_n_fold = 10 #iterations for cross-validation
	def_nu = 0.5
	def_outf = 'decision.model'

	def_w_m = [None, def_c_min, def_c_high, def_Gamma_min, def_Gamma_high, def_c_num, def_g_num, def_n_fold, def_nu, def_outf]
	return def_w_m

def get_epsilon():
    return 0.1

###########################Defaults: testset_creation.py################

def get_testset_cr_defaults():
	default_category = '1'
	default_minident = 50
	default_maxident = 98
	default_k = 0.2
	default_sampleSize = 100
	default_simulateSampleSize = 1
	default_num_processes = multiprocessing.cpu_count()
	default_mute = False
	default_window_size = 120
	default_step = 40

	def_test_set_cr = [default_category, default_minident, default_maxident,True,  default_k, default_window_size, default_num_processes, default_mute, default_step]
	return def_test_set_cr
