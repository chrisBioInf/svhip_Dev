#This module exists for simpler access to configuration files for scaling
# - i.e. the txt files containg the range of z-score, SCI, Shannon-Entropy
#for a model needed for normalization to interval [-1, 1].
#Figured it might be useful to be able to edit/add values manually.
#
#scale.config has exactly 1 line in format 'z-score low : z_score high : 
#SCI low: SCI high : Shannon low : Shannon high'

def get_z_range(filename):
    try:
        with open(filename,'r') as scal:
            for line in scal.readlines():
                return [float(line.split(':')[0]), float(line.split(':')[1])]
    except Exception as e:
        raise e

def get_SCI_range(filename):
    try:
        with open(filename,'r') as scal:
            for line in scal.readlines():
                return [float(line.split(':')[2]), float(line.split(':')[3])]
    except Exception as e:
        raise e
        
def get_shannon_range(filename):
    try:
        with open(filename,'r') as scal:
            for line in scal.readlines():
                return [float(line.split(':')[4]), float(line.split(':')[5])]
    except Exception as e:
        raise e

def write_scale(data_vector, filename):
    z_values = []
    SCI_values = []
    sh_values = []
    with open(filename, 'a') as scal:
        for i in range(0, len(data_vector)):
            if None not in data_vector[i]:
                z_values.append(data_vector[i][0])
                SCI_values.append(data_vector[i][1])
                sh_values.append(data_vector[i][2])
            else:
                print('Faulty alignment data tuple detected. Skipping...')
        
        if len(z_values) == 0:
            z_values = [0]
        if len(SCI_values) == 0:
            SCI_values = [0]
        if len(sh_values) == 0:
            sh_values = [0]

        min_z, max_z = min(z_values), max(z_values)
        min_SCI, max_SCI = min(SCI_values), max(SCI_values)
        min_sh, max_sh = min(sh_values), max(sh_values)
        scal.write(str(min_z)+' : '+str(max_z)+' : '+str(min_SCI)+' : '+str(max_SCI)+' : '+str(min_sh)+' : '+str(max_sh) + '\n')
    return [min_z, max_z, min_SCI, max_SCI, min_sh, max_sh]
