import sys
import os

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
# sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'data_gen/')))
# sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'scaling/python/')))
# sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'write_m/')))
# sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'currysoup/')))
sys.path.append(os.path.abspath(os.path.join(THIS_FOLDER, 'logger/')))

import svhip.data_gen.generate_data as generate_data
import svhip.write_m.write_model as write_model
# import scaling.python.scaling_config as scaling_config
from currysoup.currysoup import *
import logger


class Svhip:

    inputfile = ""
    outputfile = ""
    argx = []
    func_argx = []
    logfile = None

    def __init__(self, argx=[]):
        print("Welcome.")
        self.argx = argx

        if '-ldcfg' in self.argx:
            try:
                self.func_argx = load_conf(self.argx[self.argx.index('-ldcfg') +1])
            except Exception as e:
                raise e
        if help_needed(self.argx, self.func_argx):
            self.help()

        if ('-h' not in self.argx and '--man' not in self.argx) and ('-h' not in self.func_argx and '--man' not in self.func_argx):
            self.function_log = logger.log_handler(self.argx, self.func_argx)
        else:
            self.function_log = None
        self.inputfile = designate_inputfile(self.argx, self.func_argx, self.function_log)
        self.outfile = designate_outputfile(self.argx, self.func_argx, self.inputfile, self.function_log)

        if '-svcfg' in argx:
            self.save_argx()

    def create_data_set(self):
        if len(self.func_argx) == 0:
            self.func_argx = self.argx
        generate_data.testset_cr(self.func_argx, self.inputfile, self.outfile, self.function_log)

    def train_model(self):
        if len(self.func_argx) == 0:
            self.func_argx = self.argx
        write_model.write_m(self.func_argx, self.inputfile, self.outfile, self.function_log)

    def verify_model(self):
        pass

    def get_input_path(self):
        return self.inputfile

    def get_output_path(self):
        return self.outputfile

    def get_argx(self):
        return self.argx

    def get_func_argx(self):
        return self.func_argx

    def save_argx(self):
        save_configuration(self.argx, self.func_argx, self.inputfile)

    def help(self):
        with open(os.path.join(THIS_FOLDER, "help/general.help"), "r") as fh:
            print(fh.read())
        return None

    def print_log(self):
        try:
            f = open(self.logfile.get_log_name(), 'r')
            print(f)
            f.close()
        except Exception as e:
            print("WARNING: Log file could not be opened. ")

            if self.logfile == "":
                print("---> Logfile seems not to be initialized. Please contact developer.")
            else:
                print("---> Please check file in location: " + str(self.logfile.get_log_name()))
            raise e
