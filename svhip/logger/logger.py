import os
import sys
import logging
import time
import traceback

LOGGER_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.abspath(os.path.join(LOGGER_DIR, 'logs/'))
TMP_LOG_DIR = '' #current writing directory - defined by output file location
PAR_DIR = os.path.abspath(os.path.join(LOGGER_DIR, os.pardir))

class log_handler:
    def __init__(self, argx, func_argx = None):
        self.log_name = self.init_logfile(argx, func_argx)
        self.tmp_log_name = self.designate_secondary_log_dir(argx)    
    def write_warning(self, string):
        logging.basicConfig(filename=self.get_log_name(),level=logging.DEBUG)
        logging.warning(string)
        logging.basicConfig(filename=self.get_tmp_log_name(),level=logging.DEBUG)
        logging.warning(string)
    
    def write_log(self, string = None):
        if string is not None:
            logging.basicConfig(filename=self.get_log_name(),level=logging.DEBUG)
            logging.debug(string)
            logging.basicConfig(filename=self.get_tmp_log_name(),level=logging.DEBUG)
            logging.debug(string)
        
        with open(self.get_log_name(), 'a') as lg:
            lg.write(string + '\n')
        self.begin_traceback(self.get_log_name())
        
        with open(self.get_tmp_log_name(), 'a') as lg:
            lg.write(string + '\n')
        self.begin_traceback(self.get_tmp_log_name())
        
    def begin_traceback(self, logfile):
        save_stdout = sys.stdout
        
        with open(logfile, 'a') as lg:
            lg.write('MOST RECENT TRACEBACK: \n')
            sys.stdout = lg
            rep_lg = traceback.format_stack()
            print(' '.join(rep_lg))
        sys.stdout = save_stdout
    
    def init_logfile(self, argx, func_argx):
        '''
        Initialize a simple logfile, that is passed down with remaining arguments.
        It will be filled with system tracebacks, warnings and debug messages.
        Can optionally be deleted on a succesfull run.
        '''
        t = time.localtime()
        now = time.strftime("%H:%M:%S", t)
        log_name = LOG_DIR + '/log_' + str(now)
        argx = argx[1:]
        
        with open(log_name, 'w+') as lg:
            lg.write('System time: ' + str(now) +'\n')
            lg.write('Argument command line: ' + ' '.join(argx) +'\n') 
            if func_argx is not None and len(func_argx)> 0:
                lg.write('Argument loaded from config file: ' + ' '.join(func_argx) +'\n')
            lg.write('\n')
        
        return log_name
    
    def designate_secondary_log_dir(self, argx):
        if '-o' in argx and (argx.index('-o')+1 < len(argx)):
            out_arg = argx[argx.index('-o')+1]
            if os.path.isdir(out_arg):
                return os.path.abspath(os.path.join(out_arg, 'svhip_logs/'))
            else:
                return os.path.abspath(os.path.join(os.getcwd(), 'svhip_logs/'))
        else:
            return os.path.abspath(os.path.join(os.getcwd(), 'svhip_logs/'))
        
    
    def get_log_name(self):
        return self.log_name
    
    def get_tmp_log_name(self):
        return self.tmp_log_name
