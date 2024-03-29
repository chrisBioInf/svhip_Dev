
import os

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))

from .svhip_object import Svhip
from .currysoup.currysoup import *

def main():
    '''
    First Initialization.
    TODO self: The situation of the func_argx parameter is kinda awkward right 
    now. Maybe later check restructuring of main argx at beginning of first
    call.
    '''
    svhip = Svhip(argx=sys.argv)
 
    '''
    Interpret remaining command line:
    '''
    if len(svhip.argx) <= 1:
        show_general_help()
    
    if '-data_gen' in svhip.argx or '-data_gen' in svhip.func_argx:
        svhip.create_data_set()
    
    if '-train' in svhip.argx or '-train' in svhip.func_argx:
        svhip.train_model()
    
    elif '-auto' in svhip.argx or '-auto' in svhip.func_argx:
        svhip.full_run()
    
    print('Program succesfully finished. Now closing.')
    return 1


###########################Misc.########################################

def show_general_help():
    with open(os.path.join(THIS_FOLDER, "help/general.help"), "r") as fh:
        print(fh.read())
    return None

if __name__ == "__main__":
    main()
