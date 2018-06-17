#!/usr/bin/env python

import argparse
from typing import List, Any, Union

import numpy as np
from scipy.constants import alpha
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import SGDRegressor

import PeptideScores as Ps


def evaluate_result(result: List[int], targets: np.ndarray, target_labels: List[str]):
    tp = fp = tn = fn = 0.0
    for i in range(len(result)):
        if result[i] == targets[i]:
            if result[i] == 1:
                tp += 1
            else:
                tn += 1
        else:
            print("Wrong prediction for {}: true class {}, predicted class {}".format(
                target_labels[i], targets[i], result[i])
            )
            if result[i] == 1:
                fp += 1
            else:
                fn += 1

    acc = (tp + tn) / (tp + tn + fp + fn)
    sensitivity = tp / (tp + fn)
    specificity = tn / (fp + tn)

    # Output result
    print(
        "Result:\nTP: {}\tTN: {}\nFP: {}\tFN: {}\nAccuracy: {:.3%}\tSensitivity: {:.3%}\tSpecifity: {:.3%}".format(
            tp, tn, fp, fn, acc, sensitivity, specificity
        ))
    pass


def sum_2d_list_columns(list_2d: List[List[int]]) -> List[int]:
    b: List[int] = []
    for i in range(len(list_2d[0])):
        b.append(0)
        for j in range(len(list_2d)):
            b[-1] = b[-1] + 1 if (list_2d[j][i] == 1) else -1
        b[-1] = 1 if (b[-1] / len(list_2d) >= 0.5) else 0
    return b


def main():
    print("Parsing arguments...")
    # Parse argumente
    parser = argparse.ArgumentParser(description="Random Epitope Binding Predictor")
    parser.add_argument('--train', metavar='training-file', help='path to the file to train from', required=False)
    parser.add_argument('--input', metavar='in-file', help='path to the input file', required=False)
    parser.add_argument('--output', metavar='out-file', help='path to save the output file to', required=False)
    args = parser.parse_args()
    # Default path of input file
    if args.train is None:
        args.train = "data\\project_training.txt"

    # Parse Training data in to PeptidScoring objekt
    training_set = Ps.PeptideScores()
    print("Parsing training file from ", args.train)
    try:
        file = open(args.train, "r")
        training_set.parse(file.readlines())
        file.close()
    except FileNotFoundError:
        print("Could not find training file")
        return 1
    except EOFError:
        print("Unexpected error")
        return 2

    data_set_number = 3
    training_set_observations = np.array([])
    training_set_targets = np.array([])
    testing_set_observations = np.array([])
    testing_set_targets = np.array([])
    testing_set_observation_peptides = []

    # Parse input file
    if args.input is not None:
        testing_set = Ps.PeptideScores()
        print("Parsing training file from ", args.input)
        try:
            file = open(args.input, "r")
            testing_set.parse(file.readlines())
        except FileNotFoundError:
            print("Could not find input file")

        training_set_observations = training_set.get_data(data_set_number)
        training_set_targets = training_set.get_targets()
        testing_set_observations = testing_set.get_data(data_set_number)
        testing_set_targets = testing_set.get_targets()
        testing_set_observation_peptides = testing_set.get_peptides()
    else:
        # Split the test and training data when no input file is given
        split = int(len(training_set.observations) * 1.0 / 10)
        split_train_sets = training_set.shuffle_and_split([0, split, -1], shuffle=True)

        training_set_observations = split_train_sets[1].get_data(data_set_number)
        training_set_targets = split_train_sets[1].get_targets()
        testing_set_observations = split_train_sets[0].get_data(data_set_number)
        testing_set_targets = split_train_sets[0].get_targets()
        testing_set_observation_peptides = split_train_sets[0].get_peptides()

    print("Using Data:\n",
          training_set_observations, "\n",
          training_set_targets, "\n",
          testing_set_observation_peptides, "\n",
          testing_set_observations, "\n",
          testing_set_targets)
    # Scale the training data and test data to the training data
    standard_scaler = StandardScaler()
    standard_scaler.fit(training_set_observations)
    training_set_observations = standard_scaler.transform(training_set_observations)
    testing_set_observations = standard_scaler.transform(testing_set_observations)

    print(len(training_set_observations), " elements in the training set,", len(testing_set_observations),
          "in the testing set")
    # Train the model
    classifiers = [
        RandomForestClassifier(bootstrap=False, max_depth=10, max_features='sqrt', n_estimators=80),
        MLPClassifier(alpha=0.1, hidden_layer_sizes=(500, 250, 125, 25), activation='logistic', max_iter=500),
        MLPClassifier(alpha=0.001, hidden_layer_sizes=(1000, 700, 350, 175, 25), activation='relu', max_iter=500)
    ]

    # Train the model
    result: List[List[int]] = []
    for ci in range(len(classifiers)):
        print("Training and evaluating Classifier", ci)
        # train classifier
        classifiers[ci].fit(training_set_observations, training_set_targets)

        # cross validate
        scores = cross_val_score(classifiers[ci], testing_set_observations, testing_set_targets, cv=2)

        print(scores)
        print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

        # classifier.predict()
        result.append(classifiers[ci].predict(testing_set_observations))

        # evaluation
        evaluate_result(result[ci], testing_set_targets, testing_set_observation_peptides)

    print("Evaluation Total:")
    evaluate_result(sum_2d_list_columns(result), testing_set_targets, testing_set_observation_peptides)
    return 0


if __name__ == '__main__': main()
