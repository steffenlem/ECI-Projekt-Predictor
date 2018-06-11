#!/usr/bin/env python

import argparse
from typing import List, Any
from sklearn.neighbors import KNeighborsClassifier
from sklearn.utils import Bunch
from sklearn import tree
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import pandas
import itertools
import aminoAcidEncoding as encode



class PeptideScores(object):
    feature_names = []
    observations = []
    target_names = []
    target = []

    def parse(self, string : List[str]):
        # Nimm die erste Spalte als tabellen header
        # 1. Peptid sequenz, 2. IC50 score, 3. targets
        header_tab_sep = string.pop(0).replace("\n","").split("\t")
        # speicher in die Peptid klasse
        self.add_names([header_tab_sep[0], header_tab_sep[1]], header_tab_sep[2])
        # lese jede line und splitte sie wie oben
        for line in string:
            line_tab_sep = line.replace("\n","").split("\t")
            # Nimm die Spalten als datenpunkt in der peptide class auf
            # 1. Peptid sequenz, 2. IC50 score, 3. targets
            self.add_observation([line_tab_sep[0],line_tab_sep[1]], line_tab_sep[2])
        pass

    def get_target_name(self):
        return self.target_names

    def add_names(self, feature_names, target_names):
        self.feature_names = feature_names
        self.target_names  = target_names
        pass

    def add_observation(self, observation, target):
        self.observations.append(observation)
        self.target.append(target)
        pass

    def get_bunch(self):
        return Bunch(data          = self.observations,
                     target        = self.target,
                     target_names  = self.target_names,
                     feature_names = self.feature_names)

#Determines a random test/train set for a given ratio
#Parameters
#   data            array of data
#   classificaton   array of claffifications
#   test_ratio
#
#Returns
#

def split_train_test(data, classification, test_ratio):
    shuffled_indices = np.random.permutation(len(data))
    test_set_size = int(len(data) * test_ratio)
    test_indices = shuffled_indices[:test_set_size]
    train_indices = shuffled_indices[test_set_size:]
    return [data[i] for i in train_indices], [classification[i] for i in train_indices], [data[i] for i in test_indices], [classification[i] for i in test_indices]





def main():
    print("Parsing arguments...")
    # Parse argumente
    parser = argparse.ArgumentParser(description="Random Epitope Binding Predictor")
    parser.add_argument('--train', metavar='training-file', help='path to the file to train from', required=False)
    parser.add_argument('--input', metavar='in-file', help='path to the input file', required=False)
    parser.add_argument('--output', metavar='out-file', help='path to save the output file to', required=False)
    args = parser.parse_args()
    # Verwende default path der input file
    if (args.train == None):
        args.train = "data\\project_training.txt"
    # Verwende default path der input file
    if (args.input == None):
        args.input = "data\\test_input.txt"
    # Verwende default path der input file
    if (args.output == None):
        args.output = "out\\output.txt"

    # lade instanzen der benötigten klassen
    pepscore = PeptideScores()
    classifier = KNeighborsClassifier(n_neighbors=1)
    #parse die text file
    print("Parsing file from ", args.train)
    try:
        file = open(args.train, "r")
        pepscore.parse(file.readlines())
        file.close()
    except FileNotFoundError:
        print("Could not find File")
    except:
        print("Unexpected error")

    print("Found:")
    print("Feature names: ", pepscore.feature_names)
    print("Target names: ", pepscore.target_names)
    print("Observations: ", pepscore.observations)
    print("Targets: ", pepscore.target)

    #remove ic50 and split sequence
    splitted_training_data = []
    for x in pepscore.observations:
        splitted_training_data.append(list(x[0]))

    #transform to 20 bit endcoding
    encoded_20bit_data = encode.twentybit(splitted_training_data)

    #transform to blomap
    encoded_blomap_data = encode.blomapCode(splitted_training_data)


    #split train/test
    train_set, train_set_target, test_set, test_set_target = split_train_test(encoded_blomap_data, pepscore.target,0.4)

    #print(train_set)
    #print(train_set_target)

    #transfrom data to numpy array
    train_set = np.array(train_set)
    train_set_target = np.array(train_set_target)
    test_set = np.array(test_set)
    test_set_target = np.array(test_set_target)


    # Train the model
    #clf = tree.DecisionTreeClassifier()
    #clf = MLPClassifier(solver='lbfgs', alpha=0.001,hidden_layer_sizes = (180, 100), random_state= 1, activation='relu')
    clf = RandomForestClassifier(bootstrap= False, max_depth= 10,max_features= 'sqrt',n_estimators = 80)
    scores = cross_val_score(clf, np.array(encoded_blomap_data), np.array(pepscore.target),cv=5)

    print(scores)
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))



    ##Grid Search
    # param_grid = [
    #     {'solver': ['lbfgs'], 'alpha': [0.01, 0.001, 0.0001, 0.00001], 'activation': ['relu','logistic', 'tanh', 'identity'], 'random_state': [1], 'hidden_layer_sizes':[20,25,30,35,40]},
    #     #{'solver': ['lbfgs'], 'alpha': [0.01, 0.001, 0.0001, 0.00001], 'activation': ['relu'], 'random_state': [1], 'hidden_layer_sizes':[x for x in itertools.product((40,60,120,100,140,80),repeat=2)]}
    #     #{'bootstrap': [True, False], 'max_depth': [10, 30, 50, 70, None],'max_features': ['auto', 'sqrt'],'n_estimators': [10,20,40,80,160,200,300]}
    #     #,'min_samples_leaf': [1, 2, 4],'min_samples_split': [2, 5, 10]
    #     #{'solver': ['adam', 'sgd'], 'alpha': [0.001, 0.0001, 0.00001, 0.000001],'learning_rate': ["constant", "invscaling", "adaptive"], max_iter=1000,}
    # ]
    #
    #
    #
    # grid_search = GridSearchCV(MLPClassifier(), param_grid, cv=5, scoring='accuracy')
    #
    # grid_search.fit(np.array(encoded_20bit_data), np.array(pepscore.target))
    #
    # print(" ")
    # print("Best Parameter:")
    # print(grid_search.best_params_)
    # cvres = grid_search.cv_results_
    #
    # print("")
    # print("Top 10 combinations (unordered):")
    # for mean_test_score, std_test_score, rank_test_score, params in zip(cvres["mean_test_score"], cvres["std_test_score"], cvres["rank_test_score"],cvres["params"]):
    #    if(rank_test_score in range(1,10)):
    #        print(mean_test_score, std_test_score, params)


    #Train/Test Split Part auskommentiert
    """
    clf = clf.fit(train_set, train_set_target)
    
    #classifier.predict()
    result = clf.predict(test_set)

    #evaluation
    tp = fp = tn = fn = 0
    for i in range(len(result)):
        if ((result[i] == test_set_target[i]) and (test_set_target[i] == '1')):
            tp += 1
        elif ((result[i] != test_set_target[i]) and (test_set_target[i] == '1')):
            fn += 1
        elif ((result[i] == test_set_target[i]) and (test_set_target[i] == '0')):
            tn += 1
        elif ((result[i] != test_set_target[i]) and (test_set_target[i] == '0')):
            fp += 1






    acc = (tp + tn) / (tp + tn + fp + fn)
    sensitivity = tp / (tp + fn)
    specificity = tn / (fp + tn)



    #Output result
    print()
    print("TP: " + str(tp))
    print("TN: " + str(tn))
    print("FP: " + str(fp))
    print("FN: " + str(fn))
    print("Accuracy: " + str(acc))
    """
    return 0










if __name__ == '__main__': main()