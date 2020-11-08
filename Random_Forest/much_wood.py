import numpy as np
import sklearn.ensemble as sk
from sklearn.datasets import make_classification
import matplotlib.pyplot as plt
import sklearn.metrics as m
import pickle 

def create_forest(X, Y):
    frst = sk.RandomForestClassifier(max_depth=10, random_state = 0)
    frst.fit(X, Y)
    return frst

def save_forest(frst, path):
    pickle.dump(frst, open(path, 'wb'))
    
def load_forest(path):
    return pickle.load(open(path, 'rb'))

def predict_instance(frst, X):
    try:
        return frst.predict_proba(X), frst.predict(X)
    except ValueError:
        return frst.predict_proba([X])
    print("Could not obtain valid classification. Aborting...")
    return None

def woodchopper(filename):
    X, Y, s, v = [], [], [], []
    
    with open(filename, 'r') as f:
        for line in f.readlines():
            v = []
            s = line.split(' ')
            Y.append(int(s[0]))
            v.append(float(s[1].split(':')[1]))
            v.append(float(s[2].split(':')[1]))
            v.append(float(s[3].split(':')[1]))
            X.append(v)
    
    return np.asarray(X), np.asarray(Y)

def return_metric(y_true, y_score, y_predict):
    AUC = round(m.roc_auc_score(y_true, y_score), 5)
    f1 = round(m.f1_score(y_true, y_predict), 5)
    return AUC, f1
    
def draw_roc(y_true, y_score, y_predict):
    fpr, tpr, threshold = m.roc_curve(y_true, y_score)
    del(threshold)
    AUC, f1 = return_metric(y_true, y_score, y_predict)
    fig, ax = plt.subplots()
    plt.plot(fpr, tpr, 'g')
    plt.xlabel("false positive rate")
    plt.ylabel("true positive rate")
    plt.ylim([0.0, 1.1])
    plt.xlim([-0.1, 1.0])
    plt.plot([0,1], [0,1], 'k--')
    plt.text(0.5, 0.25, "Area under curve: " + str(AUC) + " \n F1_score: " + str(f1))
    #plt.show()
    plt.savefig("forest_vs_oldset.pdf")
    
def main(trainset, testset, name):
    X_train, Y_train = woodchopper(trainset)
    X_test, y_true = woodchopper(testset)
    y_score = []
    
    frst = create_forest(X_train, Y_train)
    save_forest(frst, name)
    
    print("Now predicting...")
    frst2 = load_forest(name)
    Y_out, y_predict = predict_instance(frst2, X_test)
    for y in Y_out:
        y_score.append(y[1])
    
    draw_roc(y_true, np.asarray(y_score), np.asarray(y_predict))

main("trainingsets/sample_1.dat", "testsets/oldset.data", "example_Forest.forest")

''' ####################################################
Example code - generates minimal data set and forest 

X, Y = make_classification(n_samples = 100, n_features = 3, 
                           n_informative =2,n_redundant = 0,
                           random_state = 0, shuffle = False)

print(X)
print(Y)
frst = sk.RandomForestClassifier(max_depth=2, random_state = 0)
frst.fit(X, Y)
print(frst.predict([[10,-8,-3], [3,4,4]])[1])
'''