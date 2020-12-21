#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 09:24:23 2020

@author: christopher
"""

import numpy as np
from sklearn import svm
from sklearn.linear_model import LogisticRegression
import sklearn.decomposition as decompose
import sklearn.ensemble as sk
from sklearn.datasets import make_classification
import matplotlib.pyplot as plt
import sklearn.metrics as metrics

def initialize_classifiers():
    
    svm_classifier = svm.SVC(C = 1.0, 
                             kernel = 'rbf', 
                             probability= True )
    
    forest_classifier = sk.RandomForestClassifier(n_estimators = 500, 
                                     max_leaf_nodes=16, 
                                     random_state = 0,
                                     min_samples_split= int(0.01 * 10000),
                                     min_samples_leaf = int(0.01 * 10000),
                                     max_depth=None
                                     )
    
    logistic_regressor = LogisticRegression()
    
    return svm_classifier, forest_classifier, logistic_regressor

def build_ensemble(X, Y):
    
    svm_clf, forest_clf, logistic_clf = initialize_classifiers()
    
    majority = sk.VotingClassifier(
        estimators = [('SVM', svm_clf), ('Random Forest', forest_clf), ('Logistic', logistic_clf)],
        voting = 'hard')
    
    majority.fit(X, Y)
    
    return majority


class feature_vector:
    SCI = 0
    z_score = 0
    entropy = 0
    category = 0
    
    def scale_back(self, x, mn, mx):
        return ( (x - (-1)) / (1 - (-1)) ) * (mx - mn ) + mn   
    
    def __init__(self, feature_string):
        a, s, z, e = feature_string.split(' ')
        self.SCI = self.scale_back(float(s.split(':')[1]), 0.0, 1.29)
        self.z_score = self.scale_back(float(z.split(':')[1]), -8.15, 2.01)
        e = e.replace('\n', '')
        self.entropy = self.scale_back(float(e.split(':')[1]), 0.0, 1.29)
        self.category = -int(a)
        
    def print_vector(self):
        print("SCI: " + str(self.SCI) + ", z-Score: " + str(self.z_score) + ", Shannon-entroy: " + str(self.entropy))
        
    def get_vector(self):
        return (self.SCI, self.z_score, self.entropy)
    
    def get_category(self):
        return self.category
    
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
    
def read_file(name):
    vector_list = []
    
    with open(name, 'r') as f:
        for line in f.readlines():
            vector_list.append(feature_vector(line))
    return vector_list

def get_color():
    for c in ['r', 'g', 'b']:
        yield c

def draw_roc(fpr, tpr, AUC, f1, i):
    color = ['r', 'g', 'b']
    c = color[i]
    print(i, c)
    
    plt.figure(figsize = (8, 6), dpi = 250)
    plt.plot(fpr, tpr, color = c)
    plt.plot([0,1], [0, 1], 'k--')
    
    plt.xlabel('false positive rate')
    plt.ylabel('true positive rate')
    plt.text(0.5, 0.25, "Area under curve: " + str(AUC) + " \n F1 score: " + str(f1))
    plt.xlim((-0.1, 1))
    plt.ylim((0, 1.1))
    
    plt.savefig('clf_' + str(i) + '_roc.pdf')
    
    plt.show()

def main():
    X, Y = woodchopper('testsets/sample_1.dat')
    X_test, Y_test = woodchopper('testsets/oldset.data')
    
    triumvirate = build_ensemble(X, Y)
    
    #Y_pred = triumvirate.predict(X_test)
    
    #print(metrics.accuracy_score(Y_test, Y_pred))
    
    clf1, clf2, clf3 = initialize_classifiers()
    
    i = 0
    
    for clf in (clf1, clf2, clf3, triumvirate):
        clf.fit(X, Y)
        Y_pred = clf.predict(X_test)
        
        if clf != triumvirate:
            y_score = clf.predict_proba(X_test)
            y_score = y_score.T
            y_score[[0, 1]] = y_score[[1, 0]]
            
            fpr, tpr, threshold = metrics.roc_curve(Y_test, y_score[0])
            AUC = round(metrics.roc_auc_score(Y_test, y_score[0]), 5)
            f1 = round(metrics.f1_score(Y_test, Y_pred), 5)
            print(clf.__class__.__name__, 
              metrics.accuracy_score(Y_test, Y_pred),
              "\n Area under curve: " + str(AUC),
              "\n F1_score: " + str(f1))
            draw_roc(fpr, tpr, AUC, f1, i)
            i = i + 1
    
        else:
            f1 = round(metrics.f1_score(Y_test, Y_pred), 5)
            print(clf.__class__.__name__, 
              metrics.accuracy_score(Y_test, Y_pred),
              "\n F1_score: " + str(f1))
        
main()
    