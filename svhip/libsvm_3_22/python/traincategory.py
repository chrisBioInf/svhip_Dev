import sys
#sys.path.append('/libsvm-3.22/python/') ...get this to work, dammit.
from svmutil import *
from svm import *
from svm import __all__ as svm_all
import matplotlib.pyplot as plt
from sklearn import metrics
import statistics


def write_model(filename):
    #optionsa = [-s, -t, -d, -g,-r,-c,-n,-p,-m,-e,-h,-b,-w]
    #optionsb = [0, 2, 3, 1/200, 0,1,0.5,0.1, 100, 0.001, 1, 0, 1]
    param = svm_parameter()
    param.gamma = 0.5
    param.cross_validation = True
    #print(param)
    #print(type(param))
    
    pr = svm_read_problem(filename)
    #print(pr)
    print(pr[0])
    print(pr[1][0][1])
    trset = svm_problem(pr[0], pr[1])
    print(trset)
    cmd = '-s {svm_type} -t {kernel_type} -c {c} -g {g} -n {nu} -b 1'.format(
            svm_type = 0,
            kernel_type = 2,
            c = 32768,
            g =  8,
            nu = 0.5)
    m = svm_train(trset, cmd)
    svm_save_model(filename.replace(".dat", ".model"), m)
    
    #print(type(pr))
    #print(type(trset))
    #m = svm_train(trset, param)
    #$svm_save_model(str(str(filename) + '_model'), m)    

#write_model('heart_scale.tr')

def predict_data(filename_problem, filename_model):
    #pr = svm_read_problem(filename_problem)
    #model = svm_load_model("decision_dinucleotide.model")
    #result = svm_predict(pr[0], pr[1], m, options='-b 1')

    model = svm_load_model(filename_model)
    count = 0
    scores = []
    y_true = []
    
    x = []
    y = []
    
    with open(filename_problem, 'r') as f:
        outlines = []
        values = {}
        for line in f.readlines():
            
            s = (line.split(' '))
            print(s)
            
            values[1] = (float(s[1].split(':')[1]))
            values[2] = (float(s[2].split(':')[1]))
            values[3] = (float(s[3].split(':')[1]))
            values = [values]
            
            count += 1
            y_true.append(str(s[0]))
            a,b,val = svm_predict([int(s[0])],values,model, options = "-b 1")
            print(a)
            print(b)
            print(val)
            #print("ClassA: " + str(val[0][0]) + ", ClassB: " + str(val[0][1]))
            scores.append(val[0][0])
            #x.append(statistics.mean([values[0][1], values[0][2], values[0][3]]))
            #y.append(val[0][0])
            
            outlines.append(line)
            outlines.append("ClassA: " + str(val[0][0]) + ", ClassB: " + str(val[0][1]) + '\n')
        
    fpr, tpr, thresholds = metrics.roc_curve(y_true, scores, pos_label=1)
    fig, ax = plt.subplots()

    #plt.title("Predictive power of original RNAz 2.0 model \n  on custom test set of 200 families \n and original testset from supplement")
    #plt.legend(loc = "lower right")
    plt.plot(fpr, tpr, 'b')
    plt.xlim([-0.1, 1.1])
    plt.ylim([-0.1, 1.1])
    plt.plot([-0.1, 1.1],[-0.1, 1.1], 'k--')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    
    area = calcAreaUnderCurve(fpr, tpr)
    f = metrics.f1_score(y_true, y_predict)
    #plt.text(0.45, 0.3, "Area under curve old testset: " + str(area))
    
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(name + ".pdf")
    
                
def calcAreaUnderCurve( xvalues, yvalues):
        area = 0.0
        for i in range(0, len(xvalues)-1):
            area = abs(area + calcTrapez(xvalues[i], xvalues[i+1], yvalues[i], yvalues[i+1]))
        return round(area, 3)
    
def calcTrapez( x1, x2, y1, y2):
        a = 0.0
        a = abs(a + ((x2-x1)*((y2+y1)/2)))
        return a
    
    
    #for j in range(1,121):
        #print(j)
        #a,b,val = svm_predict([1],values,model)
        #print(str(a) + ', ' + str(b) + ', ' + str(val))
        #newval = val[0][0]
        '''
        for i in range(1,5):
            values[0][i] = values[0][i+1]
        values[0][5] = newval
        '''
    
    #print(sorted([m, m])[0])

write_model("/scr/k70san2/christopher/grid_k70/svhip/libsvm_3_22/python/rfam_scaled.dat")
    # Data ; Model
#predict_data('rfam_scaled.dat', 'Vol1_newTry.model', 23703)


