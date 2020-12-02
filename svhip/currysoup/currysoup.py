import sys
import os

CURRY_DIR = os.path.dirname(os.path.abspath(__file__))
PAR_DIR = os.path.abspath(os.path.join(CURRY_DIR, os.pardir)) 
sys.path.append(os.path.abspath(os.path.join(PAR_DIR, 'logger/')))

'''
This module serves as an extension of the main module to drastically unclog
the (before barely readable) code used to interpret the command line. In 
service of maintanence it would be desirable to keep this code as ordered 
and modularised as possible.

Other than parsing command line in sizable bits, it is remarkably boring.

The origin of the name is a secret to everybody.
'''
##########################Functions called from main()##################

def help_needed(argx, func_argx):
    if ('-h' in argx and '-write_m' not in argx and '-data_gen' not in argx) or ('--man' in argx and '-write_m' not in argx and '-data_gen' not in argx) :
        return True
    elif ('-h' in func_argx or '--man' in func_argx) and '-write_m' not in func_argx and '-data_gen' not in func_argx:
        return True


def designate_inputfile(argx, func_argx, flog):
    if (('-i' not in argx) or (argx.index('-i')+1 >= len(argx))) and ('-h' not in argx and '-h' not in func_argx and '--man' not in argx and '--man' not in func_argx):
        lg_message = "FATAL: Input parameter is required. Use -h for more info."
        e = ValueError(lg_message)
        flog.write_log(lg_message)
        raise e
    
    elif (('-i' not in argx) or (argx.index('-i')+1 >= len(argx))):
        return None
    
    elif os.path.isdir(argx[argx.index('-i')+1]):
        return argx[argx.index('-i')+1]
    
    else:
        proc_path = os.getcwd()
        try:
            return os.path.abspath(os.path.join(proc_path, argx[argx.index('-i')+1]))
        except Exception as e:
            flog.write_log(str(e))
            raise e


def designate_outputfile(argx, func_argx, inputname, flog):
    if ('-o' in argx) and (argx.index('-o')+1 > len(argx)-1):
        lg_message = "FATAL: Output parameter is set but no argument given. View -h for more info."
        e = ValueError(lg_message)
        flog.write_log(lg_message)
        raise e
    
    elif ('-o' in argx) and (argx.index('-o')+1 < len(argx)):
        out_arg = argx[argx.index('-o')+1]
        if os.path.isdir(argx[argx.index('-o')+1]):
            return argx[argx.index('-o')+1]
        else:
            proc_path = os.getcwd()
            return os.path.join(proc_path, out_arg)
    
    #If there is no '-o' designation, use default as fallback:
    elif inputname is not None: 
        try:
            proc_path = os.getcwd()
            input_shlex = inputname.split('/')
            return os.path.join(proc_path, input_shlex[len(input_shlex)-1])
        except Exception as e:
            flog.write_log(str(e))
            raise e
    else:
        return None
    
def save_configuration(argx, func_argx, inputname):
    if len(argx) < argx.index('-svcfg')+1 or is_parameter(argx[argx.index('-svcfg')+1]):
        input_shlex = inputname.split('/')
        cfgname = input_shlex[-1]
        cfgname = cfgname.replace('.dat', '')
        save_conf(func_argx, cfgname)
    else:
        save_conf(func_argx, str(argx[argx.index('-svcfg') +1]))
        

##############Interpretation of testset creation parameters#############

def creation_soup(argx, inputfile, outfile, flog):
    parameters = [None]*9
    calls = []
    
    '''
    Parameters[3] contains the variable for automatic negative set generation...
    '''
    parameters[3] = True
    
    if '-readall' in argx:
        calls.append('-readall')
    if '-extract' in argx:
        calls.append('-extract')
    
    for a in range(0, len(argx)-1):
        if argx[a] == '-neg':
            parameters[0] = '-1'
            
        if argx[a] == '-minid': 
            try:
                parameters[1] = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -minid. Type has to be castable to Integer."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
    
        elif argx[a] == '-maxid':
            try:
                parameters[2] = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -maxid. Type has to be castable to Integer."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
        
        elif argx[a] == '--no-negative':
            parameters[3] = False
        
        elif argx[a] == '-k':
            try:
                parameters[4] = float(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -k. Type has to be castable to float."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
                
        elif argx[a] == '-w':
            try:
                parameters[5] = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign parameter -w (windowsize). Value has top be castable to Int."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
        
        elif argx[a] == '-nproc':
            try:
                parameters[6] = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -nproc. Type has to be castable Integer."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
                
        elif argx[a] == '-mute':
            parameters[7] = True
            
        elif argx[a] == ' -s ':
            try:
                parameters[8] = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign parameter -s (stepsize). Has to be castable to Int."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
            
        elif argx[a] == '-o':
            if a+1 < len(argx) and isinstance(argx[a+1], str) and outfile == None: 
                outfile = str(argx[a+1])
            else:
                raise TypeError("Inavlid output filename or outfile already given.")
        elif argx[a] == '-i':
            if a+1 < len(argx) and isinstance(argx[a+1], str) and inputfile == None: 
                inputfile = str(argx[a+1])
            else:
                raise TypeError("Inavlid input filename or inputfile already given.")
        else:
            continue
    
    return parameters, calls


##############Interpretation of write model parameters##################

def write_soup(argx, inputfile, outfile, flog):
    '''
    Initialize parameters.
    '''
    options = [None]*10
    nproc = None
    recursive_grid = False
    mute, epsilon = False, 0.0
    if outfile != None:
        options[9] = outfile

    '''
    Parse command line arguments applicable here:
    '''
    for a in range(0, len(argx)-1):
        if argx[a] == '-c_low':
            try:
                options[1] = float(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -c_low. Type has to be (castable to) float."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
                
        elif argx[a] == '-c_high':
            try:
                options[2] = float(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -c_high. Type has to be (castable to) float."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
                
        elif argx[a] == '-g_low':
            try:
                options[3] = float(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -g_low. Type has to be (castable to) float."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
  
        elif argx[a] == '-g_high':
            try:
                options[4] = float(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -g_high. Type has to be (castable to) float."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e

        elif argx[a] == '-num_c':
            try:
                options[5] = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -num_c. Type has to be (castable to) Integer."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e

        elif argx[a] == '-num_g':
            try:
                options[6] = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -num_g. Type has to be (castable to) Integer."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e

        elif argx[a] == '-n_fold':
            print(str(argx[a+1]))
            try:
                options[7] = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -n_fold. Type has to be (castable to) Integer."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e

        elif argx[a] == '-nu':
            try:
                options[8] = float(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -nu. Type has to be (castable to) float."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
                
        elif argx[a] == '-nproc':
            try:
                nproc = int(argx[a+1])
            except Exception:
                lg_message = "Could not assign input parameter -nproc. Type has to be (castable to) Integer."
                e = TypeError(lg_message)
                flog.write_log(lg_message)
                raise e
                
        elif argx[a] == '-e':
            try:
                epsilon = float(argx[a+1])
            except Exception as e:
                print("Could not assign input parameter -e (epsilon). Type has to be (castable to) float.")
                raise e
                
        elif argx[a] == '-r':
            recursive_grid = True

        elif argx[a] == '-mute':
            mute = True
            
        elif argx[a] == '-o':
            if outfile == None:
                try:
                    options[9] = str(argx[a+1])
                except Exception as e:
                    raise e
        elif argx[a] == '-i':
            if inputfile == None:
                try:
                    inputfile = str(argx[a+1])
                except Exception as e:
                    raise e
    
    return options, nproc, mute, epsilon, recursive_grid

        
############################Helper######################################        

def save_conf(func_argx, name):
    #THIS FEATURE IS STILL EXPERIMENTAL - BuT IT ALLOWS TO SAVE A CMD LINE CONFIGURATION FOR LATER REUSE, TO AVOID HAVING TO INPUT A LONG STRING OF PARAMETERS AGAIN

    argx = func_argx
    if '.cfg' not in name:
        name = name + '.cfg'
    cfgpath = os.path.abspath(os.path.join(PAR_DIR, 'configs/' +name))
    
    with open(cfgpath, 'w+') as cfg:
        cmdline = ''
        arg_included = False
        for a in range(0, len(argx)):
            if arg_included == False:
                if argx[a] in ('-h', '-mute', '-sim', '-write_m', '-testset', '-neg'):
                    cmdline = cmdline + str(argx[a]) + ';'
                elif argx[a] in ('-c_low', '-c_high', '-g_low', '-g_high', 
                                    '-num_c', '-num_g', '-n_fold', '-nu', 
                                    '-nproc', '-smpl', '-simsmpl', '-maxid',
                                    '-minid'):
                    cmdline = cmdline + argx[a]+';'+argx[a+1]+';'
                    arg_included = True
            else:
                arg_included = False
                continue
        
        cfg.write(str(cmdline))
    print('Parameters written to: ' + str(cfgpath))    
    return None
    
def load_conf(name):
    cfgpath = os.path.abspath(os.path.join(PAR_DIR, 'configs/' +name))
    
    try: 
        with open(cfgpath, 'r') as cfg:
            cmdline = cfg.readlines()[0]
            argx = str(cmdline).split(';')[:-1]
        print(argx)
        return argx
    except Exception as e:
        raise e
        
def is_parameter(arg):
    if arg in ('-c_low', '-c_high', '-g_low', '-g_high', 
                    '-num_c', '-num_g', '-n_fold', '-nu', 
                    '-nproc', '-smpl', '-simsmpl', '-maxid',
                    '-minid', '-h', '-mute', '-sim', '-write_m', 
                    '-testset', '-neg'):
        return True
    else:
        return False

