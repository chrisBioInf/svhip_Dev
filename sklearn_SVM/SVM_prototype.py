#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 11 13:17:44 2020

@author: christopher
"""
'''
NOTE: 
    Apparentlly, SVM implementation is internally even based on LibSVM?
    This 'is' slightly confusing and will have to research more. 
'''


from sklearn import svm 
import pickle
import sklearn.metrics as metrics
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score
from decimal import Decimal
import multiprocessing
from time import time

def parse_data_file(filename):
    Y, features = [], []
    
    with open(filename, 'r') as f:
        for line in f.readlines():
            s = line.split(' ')
            Y.append(int(s[0]))
            
            s[3].replace('\n', '')
            v = [float(s[1].split(':')[1]), float(s[2].split(':')[1]), float(s[3].split(':')[1])]
            features.append(v)
        
    return np.asarray(features), np.asarray(Y)

def call_classifier(m, X):
    return m.predict_proba(X), m.predict(X)

def train_classifier(X, Y, cost = 1.0, g = 0):
    if g <= 0 :
        m = svm.SVC(C = cost, kernel = 'rbf', probability= True )
    else:
        m = svm.SVC(C = cost, kernel = 'rbf', gamma = g, probability= True )
    m.fit(X, Y)
    
    return m
    
def save_classifier(m, path):
    pickle.dump(m, open(path, 'wb'))

def load_classifier(m):
    return pickle.load(open(m, 'wb'))

def return_metric(y_true, y_score, y_predict):
    AUC = round(metrics.roc_auc_score(y_true, y_score), 5)
    f1 = round(metrics.f1_score(y_true, y_predict), 5)
    return AUC, f1

def draw_roc(y_true, y_score, y_predict, name):
    fpr, tpr, threshold = metrics.roc_curve(y_true, y_score)
    del(threshold)
    AUC, f1 = return_metric(y_true, y_score, y_predict)
    fig, ax = plt.subplots()
    plt.plot(fpr, tpr, 'r')
    plt.xlabel("false positive rate")
    plt.ylabel("true positive rate")
    plt.ylim([0.0, 1.1])
    plt.xlim([-0.1, 1.0])
    plt.plot([0,1], [0,1], 'k--')
    plt.text(0.5, 0.25, "Area under curve: " + str(AUC) + " \n F1 score: " + str(f1))
    plt.savefig(name + '.pdf', dpi = 300)
    plt.show()

def score_weighted(scores):
    def square(x):
        return x*x
    return sum(map(square, scores ))/len(scores)

def crossvalidation(X, Y, cost, g):
    m = svm.SVC(kernel = 'rbf', C = cost, gamma = g)
    #scores = cross_val_score(m, X, Y, cv=10, scoring="f1_macro")
    scores = cross_val_score(m, X, Y, cv=10)
    return score_weighted(scores)

def unpack_crossvalidation(argx):
    X, Y, cost, g = argx
    print("Crossvalidation on pair: C=" + str(cost) + ', gamma=' + str(g))
    return crossvalidation(X, Y, cost, g)

def get_values(low, high, n):
    rel_pos = []
    d = (2**n) - 2

    for i in range(1, n +1):
        rel_pos.append((2**i )/d)

    def values(low, up, num, rel_pos):
        l = Decimal(str(low))
        u = Decimal(str(up))
        s = float(u-l)
        vals = [float((l)+Decimal(rel_pos[i]*s)) for i in range(int(num))]
        return vals

    valuelist = values(low, high, n, rel_pos)

    return np.asarray(valuelist)

def rudimentary_grid(a, b):
    points = []
    
    for c in a:
        for d in b:
            points.append((c, d))
    return points

def grid_search(X, Y, points):
    scores = []
    
    with multiprocessing.Pool(multiprocessing.cpu_count()) as mp:
        argx = []
        
        for pair in points: 
            c, g = pair
            argx.append((X, Y, c, g))
        scores = mp.map(unpack_crossvalidation , argx)
            
    return scores
    
def main():
    features, Y = parse_data_file("testsets/sample_1.dat")
    
    C = get_values(1.0, 2**15, 20)
    gamma = (0.5, 32.0, 20)
    
    start = time()
    
    points = rudimentary_grid(C, gamma)
    scores = grid_search(features, Y, points)
    
    max_index = scores.index(max(scores))
    max_pair = points[max_index]
    optimal_c, optimal_g = max_pair
    m = train_classifier(features, Y, optimal_c, optimal_g)
    
    t = time() - start
    
    problem_X, problem_Y = parse_data_file("testsets/oldset.data")
    probabilities, predicted_Y = call_classifier(m, problem_X)
    class_likelihood = [] 
    for y in probabilities:
        class_likelihood.append(y[1])
    class_likelihood = np.asarray(class_likelihood)
        
    AUC, f1 = return_metric(problem_Y, class_likelihood, predicted_Y)
    draw_roc(problem_Y, class_likelihood, predicted_Y, 'example_sklearn')
    print("Time: " + str(t) + ' s')

main()