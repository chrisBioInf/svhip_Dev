import sys
import os
import logging

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'data_gen/')))
sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'scaling/python/')))
sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'write_m/')))
sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'currysoup/')))
sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'logger/')))

import generate_data
import write_model
import scaling_config
from currysoup import *
import logger

def main():
    '''
    First Initialization.
    TODO self: The situation of the func_argx parameter is kinda awkward right 
    now. Maybe later check restructuring of main argx at beginning of first
    call.
    '''
    print("Welcome.")
    inputfile = None
    argx = sys.argv        
    func_argx = []
    '''
    Should a previous parameter configuration be loaded?
    EXPERIMENTAL FEATURE. You have been warned.
    '''
    if '-ldcfg' in argx:
        try:
            func_argx = load_conf(argx[argx.index('-ldcfg') +1])     
        except Exception as e:
            raise e
    '''
    Check for call of help function:
    --> Uses currysoup module.
    '''
    if help_needed(argx, func_argx):
        show_general_help()
        return None
    
    '''
    Initialize a logger:    
    '''
    if ('-h' not in argx and '--man' not in argx) and ('-h' not in func_argx and '--man' not in func_argx):
        function_log = logger.log_handler(argx, func_argx)
    else:
        function_log = None
            
    '''
    Designate input:
    None needed if -h call.
    '''
    inputfile = designate_inputfile(argx, func_argx, function_log)
    
    '''
    Designate result directory - either relative in package directory or 
    with absolute path via -oabs:
    '''    
    outfile = designate_outputfile(argx, func_argx, inputfile, function_log)
    
 
    '''
    Interpret remaining command line:
    '''
    if '-data_gen' in argx or '-data_gen' in func_argx:
        if len(func_argx) == 0:
            func_argx = argx[argx.index('-data_gen')+1:]
        generate_data.testset_cr(func_argx, inputfile, outfile, function_log)
    
    elif '-write_m' in argx or '-write_m' in func_argx:
        if len(func_argx) == 0:
            func_argx = argx[argx.index('-write_m')+1:]
        write_model.write_m(func_argx, inputfile, outfile, function_log)
    
    elif '-auto' in argx or '-auto' in func_argx:
        
        if len(func_argx) == 0:
            func_argx = argx[argx.index('-auto')+1:]
        
        generate_data.testset_cr(func_argx, inputfile, outfile + '.dat', function_log)
        write_model.write_m(func_argx, outfile + '.dat', outfile, function_log )
    
    '''
    Should the used parameters be saved in .configs?
    '''
    if '-svcfg' in argx:
        save_configuration(argx, func_argx, inputfile)
    
    print('Program succesfully finished. Now closing.')
    return None

###########################Misc.########################################

        
def show_general_help():
    with open(os.path.join(THIS_FOLDER, "help/general.help"), "r") as fh:
        print(fh.read())
    return None

if __name__ == "__main__":
    main()
