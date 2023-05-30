# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 11:23:25 2023

@author: Nimisha
"""

from random import random
from sklearn.metrics import accuracy_score,
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, average_precision_score
from sklearn.model_selection import RepeatedKFold
import pandas as pd
import statistics
from sklearn.svm import SVC
import collections
import numpy as np
import warnings
import copy
import time
import random
from sklearn import metrics
warnings.filterwarnings('ignore')

#change the file accordingly
df_new = pd.read_csv('TCGA_BRCA_gene_towork_new_withgenes_reduced_Basalagainstall_LASSO.csv')

df_dataset = df_new.iloc[:, :-1]
df_label = df_new.iloc[:, -1]


PV = 4

# change the C value according to subtype
svm_model = SVC(C=6)

c = 2
def generate_candidate(vector,countval):
    """
    Generates a new candidate solution based on the probability vector
    """
    value = []
    for p in vector:
        if random.random() < p:
            value.append(1) 
        else: 
            value.append(0) 
    
    if(all(v == 0 for v in value)):
        countval = countval + 1
        if (countval > 5):
            return value
        else:
            value = generate_candidate(vector,countval)
    
    return value


def generate_vector(size):
    """
    Initializes a probability vector with given size
    """
    return [0.5] * size

def fitness_function(chosen_candidate, X_train, X_test, y_train, y_test):
    X_final_test = pd.DataFrame()
    X_final_train = pd.DataFrame()
    for i in range(0,X_test.shape[1]):
        if chosen_candidate[i] == 1:
            X_final_train[X_train.columns[i]] = X_train.iloc[:,i]
            X_final_test[X_test.columns[i]] = X_test.iloc[:,i]
    

    svm_model.fit(X_final_train, y_train)
    y_pred_svm = svm_model.predict(X_final_test)
    ar = accuracy_score(y_test, y_pred_svm)
    pr = precision_score(y_test, y_pred_svm, average='weighted')
    f1 = f1_score(y_test, y_pred_svm, average='weighted')
    recall = recall_score(y_test, y_pred_svm, average='weighted')
    f = metrics.roc_auc_score(y_test, y_pred_svm)
    return f, ar, pr, f1, recall, X_final_train, X_final_test
    
def compete(s,vector, X_train_new, X_test_new, f, ac, p, f1, r):
    """
    Returns a tuple with the winner solution
    """
    winner_fitness = max(f)
    max_index = f.index(winner_fitness)
    winner_acc = ac[max_index]
    winner_pr = p[max_index]
    winner_re = r[max_index]
    winner_f1 = f1[max_index]
    winner_v  =  vector[max_index]
    X_winner  = X_train_new[max_index]
    X_test_winner  = X_test_new[max_index]
    return winner_fitness, winner_acc, winner_pr, winner_re, winner_f1, winner_v, X_winner, X_test_winner, max_index
    
    


def update_vector_local(s, vector, winner_f, population_size, set_dif, w_i):
    for i in range(len(vector)-1):
        for j in range(len(vector[0])):
            "winner vs loser"
            if s[w_i][j] != s[set_dif[i]][j]: 
             if s[w_i][j] == '1':
                vector[w_i][j] += 1.0 / float(population_size)
             else:
                 vector[w_i][j] -= 1.0 / float(population_size)
    return vector


def run(generations, X_train, X_test, y_train, y_test, population_size):
    vector = []
    
    PV_list = []
    size = X_train.shape[1]
    # this is the probability for each solution bit be 1
    for i in range(0,PV):
        vector.append(generate_vector(size))
        PV_list.append(i)
        
        
        
    flag = 0
    fitness = 0
    accuracy = []
    countval = 0
    best_old = []
    fitness_old = []
    vector_winner_old = []
    best_test_old = []
    accuracy1_old =[]
    precision1_old =[]
    recall1_old = []
    f1score1_old = []
    for i in range(generations):
        init = 0
        # print(i)
        s = []
        f = []
        ac = []
        p = []
        r = []
        f1 = []
        X_train_new = []
        X_test_new = []
        for j in range(0,PV):
            val = generate_candidate(vector[j],countval)
            # print(val)
            if(all(v == 0 for v in val)):
                final_val = copy.deepcopy(best_old)
                final_test = copy.deepcopy(best_test_old)
                final_auc =  copy.deepcopy(fitness_old)
                final_accuracy = copy.deepcopy(accuracy1_old)
                final_precision = copy.deepcopy(precision1_old)
                final_recall = copy.deepcopy(recall1_old)
                final_f1score = copy.deepcopy(f1score1_old)
                final_vector = copy.deepcopy(vector_winner_old)
                init = 1
                break
            else:
                s.append(val)
            
                fit, ar, pr, f_1, recall, X_train_val, X_test_val = fitness_function(s[j], X_train, X_test, y_train, y_test)
                f.append(fit)
                ac.append(ar)
                p.append(pr)
                f1.append(f_1)
                r.append(recall)
                X_train_new.append(X_train_val)
                X_test_new.append(X_test_val)
                
        if (init == 1):
            break
        wi = []
        # let them compete, so we can know who is the best of the pair
        winner_fitness, winner_acc, winner_pr, winner_re, winner_f1, winner_v, X_winner, X_test_winner, winner_index = compete(s, vector, X_train_new, X_test_new, f, ac, p, f1, r)
        wi.append(winner_index)
        
        set_dif = list(set(PV_list).symmetric_difference(set(wi)))
        
        
        vector = update_vector_local(s, vector, winner_fitness, population_size, set_dif, winner_index)
       
        
        for j in range(0,PV):
            temp1 = list(np.subtract(np.array(vector[winner_index]),np.array(vector[j])))
            temp2 = [i * c for i in temp1]
            temp3 = np.array(vector[j]) + np.array(temp2)
            vector[j] = temp3.tolist()
            
        
        accuracy.append(winner_fitness)
        
        
        if flag == 1:
            if winner_fitness > fitness:
                best = copy.deepcopy(X_winner)
                accuracy1 = copy.deepcopy(winner_acc)
                precision1 = copy.deepcopy(winner_pr)
                recall1 = copy.deepcopy(winner_re)
                f1score1 = copy.deepcopy(winner_f1)
                best_test = copy.deepcopy(X_test_winner)
                fitness = copy.deepcopy(winner_fitness)
                vector_winner = copy.deepcopy(vector[winner_index])
        
        else:
            best = X_winner
            accuracy1 = copy.deepcopy(winner_acc)
            precision1 = copy.deepcopy(winner_pr)
            recall1 = copy.deepcopy(winner_re)
            f1score1 = copy.deepcopy(winner_f1)
            best_test = X_test_winner
            fitness = winner_fitness
            vector_winner = copy.deepcopy(vector[winner_index])
            flag = 1 
        
        best_old = copy.deepcopy(best)
        best_test_old = copy.deepcopy(X_test_winner)
        accuracy1_old = copy.deepcopy(accuracy1)
        precision1_old = copy.deepcopy(precision1)
        recall1_old  = copy.deepcopy(recall1)
        f1score1_old  = copy.deepcopy(f1score1)
        fitness_old = copy.deepcopy(fitness)
        vector_winner_old = copy.deepcopy(vector_winner)
        if i > 5:
            res = accuracy[-5:]
            chk = all(x==res[0] for x in res)
            if chk:
            # if vector == v_old:
                final_val = copy.deepcopy(best)
                final_test = copy.deepcopy(best_test)
                final_auc =  copy.deepcopy(fitness)
                final_accuracy = copy.deepcopy(accuracy1)
                final_precision = copy.deepcopy(precision1)
                final_recall = copy.deepcopy(recall1)
                final_f1score = copy.deepcopy(f1score1)
                final_vector = copy.deepcopy(vector_winner)
                
                break
        if (i == 999):
            final_val = copy.deepcopy(best)
            final_test = copy.deepcopy(best_test)
            final_auc =  copy.deepcopy(fitness)
            final_accuracy = copy.deepcopy(accuracy1)
            final_precision = copy.deepcopy(precision1)
            final_recall = copy.deepcopy(recall1)
            final_f1score = copy.deepcopy(f1score1)
            final_vector = copy.deepcopy(vector_winner)
        
            
    return final_val, final_test, final_auc, final_accuracy, final_precision, final_recall, final_f1score, final_vector
        





if __name__ == '__main__':
    final_val = []
    final_val_test = []
    final_accuracy = []
    final_precision =[]
    final_recall = []
    final_f1score =[]
    final_auc = []
    final_vector = []
    row,col = df_dataset.shape
    begin = time.time()
    for i in range(0,50):
        print(i)
        training_indices = random.sample([i for i in range(row)], int(row*0.8))
        test_indices = list(set([i for i in range(row)])-set(training_indices))
        X_train = df_dataset.iloc[training_indices]
        X_test = df_dataset.iloc[test_indices]
        y_train = df_label [training_indices]
        y_test = df_label[test_indices]
        val, val_test, auc, accuracy, precision, recall, f1score, f_v = run(1000, X_train, X_test, y_train, y_test, 10)
     
        final_val.append(val)
        final_val_test.append(val_test)
        final_auc.append(auc)
        final_accuracy.append(accuracy)
        final_precision.append(precision)
        final_recall.append(recall)
        final_f1score.append(f1score)
        final_vector.append(f_v)
    end = time.time()
    fa = max(final_auc)
    
    fv = final_val[final_auc.index(fa)]
    fvtest = final_val_test[final_auc.index(fa)]
    fvector = final_vector[final_auc.index(fa)]
    faccuracy = final_accuracy[final_auc.index(fa)]
    fprec = final_precision[final_auc.index(fa)]
    frecall = final_recall[final_auc.index(fa)]
    f1score = final_f1score[final_auc.index(fa)]
    
    df_train_temp = df_new.iloc[fv.index]
    df_train = df_train_temp['Label']
    df_train_tot = pd.concat([fv, df_train],axis = 1)
    
    df_test_temp = df_new.iloc[fvtest.index]
    df_test = df_test_temp['Label']
    df_test_tot = pd.concat([fvtest, df_test],axis = 1)
    df_concat = pd.concat([df_train_tot, df_test_tot], axis=0)
    print(f'Time: {end - begin}')
print("best value:%s best value test: %s best fitness: %f" % (fv,fvtest, fa))
df_concat.to_csv('TCGA_BRCA_gene_towork_new_withgenes_reduced_Basalagainstall_mCGA.csv', index = False)
print(faccuracy)
print(fprec)
print(f1score)
print(frecall)
