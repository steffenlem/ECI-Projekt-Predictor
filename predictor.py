#!/usr/bin/env python

import argparse
from typing import List, Any
from sklearn.neighbors import KNeighborsClassifier
from sklearn.utils import Bunch


class PeptideScores:
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

    # Train the model
    #classifier.fit(pepscore.observations,pepscore.target)
    #classifier.predict()
    return 0

if __name__ == '__main__': main()