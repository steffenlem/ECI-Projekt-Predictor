from typing import List

import numpy as np

import Observation as Os


class PeptideScores(object):
    def __init__(self,
                 feature_names : List[str] = None,
                 observations  : List[Os.Observation] = None,
                 target_name   : str = None):
        """
        Creates a new instance of the Peptide scoring object
        :param feature_names: names of the features
        :param observations: observation datapoints
        :param target_names: name of the target
        """
        self.feature_names : List[str]            = [] if feature_names is None else feature_names
        self.observations  : List[Os.Observation] = [] if observations is None else observations
        self.target_name   : str                  = [] if target_name is None else target_name
        pass

    def parse(self, string: List[str]):
        """
        Parses a list of strings containing the data and saves each line as an Observation object
        :param string:
        """
        if string[0] == "Peptide\tIC50.nM.\tclass\n":
            # remove and separate the header
            header_tab_sep = string.pop(0).replace("\n", "").split("\t")
            # Save the header names
            self.add_names([header_tab_sep[0], header_tab_sep[1]], header_tab_sep[2])
        # read every line and split it by tab after removing the next line chars
        for line in string:
            line_tab_sep = line.replace("\n", "").split("\t")
            # Take every line as an Observation
            self.add_observation(peptide=line_tab_sep[0], ic50=float(line_tab_sep[1]), target=int(line_tab_sep[2]))
        pass

    def get_target_name(self) -> str:
        """
        Returns the target name
        :return:
        """
        return self.target_name

    def add_names(self, feature_names: List[str], target_name: str = None):
        """
        Adds the names for the features and the target to the wrapper object
        :param feature_names:
        :param target_name:
        """
        self.feature_names = feature_names
        self.target_name = target_name
        pass

    def get_peptides(self) -> List[str]:
        """
        Returns the peptide strings as a numpy array
        :return:
        """
        return [x.peptide for x in self.observations]

    def get_ic50(self) -> np.array:
        """
        Returns the ic50 values for all observations as a numpy array
        :return:
        """
        return np.array([x.ic50 for x in self.observations])

    def get_bit_code(self) -> np.array:
        """
        Returns the bit code for all observations as a numpy array
        :return:
        """
        return np.array([x.bit_code for x in self.observations])

    def get_blo_map(self) -> np.array:
        """
        Returns the blomap encoding for all observations as a numpy array
        :return:
        """
        return np.array([x.blo_map for x in self.observations])

    def get_blosum(self) -> np.array:
        """
        Returns the blosum encoding for all observations as a numpy array
        :return:
        """
        return np.array([x.blosum for x in self.observations])

    def get_mixed_set(self) -> np.array:
        """
        returns both blomap and bit encoding in one array
        :return:
        """
        # [x.blo_map.extend([1.0 if y else 0.0 for y in x.bit_code]) for x in self.observations]
        observation_mix : List[float] = []
        for i in range(len(self.observations)):
            a = [1.0 if y else 0.0 for y in self.observations[i].bit_code]
            a.extend(self.observations[i].blo_map)
            observation_mix.append(a)
        return observation_mix

    def get_targets(self) -> np.array:
        """
        Returns the targets for all the observations as a numpy array
        :return:
        """
        return np.array([]) if self.target_name is None else np.array([x.target for x in self.observations])

    def get_data(self, set_number : int) -> np.array:
        """
        Selects the target data set by number
        :param set_number:
        :return:
        """
        return {
            0:self.get_peptides(),
            1:self.get_bit_code(),
            2:self.get_blo_map(),
            3:self.get_mixed_set(),
            4:self.get_blosum()
        }.get(set_number,self.get_peptides())

    def add_observation(self,
                        peptide: str,
                        ic50: float,
                        bit_code: List[bool] = None,
                        blo_map: List[float] = None,
                        target: int = None):
        """
        Adds an observation to the observation list
        :param peptide:
        :param ic50:
        :param target:
        :param bit_code:
        :param blo_map:
        """
        observation = Os.Observation(peptide=peptide, ic50=ic50, bit_code=bit_code, blo_map=blo_map, target=target)
        self.observations.append(observation)

    def copy_observation(self, observation: Os.Observation):
        """
        Copies an already existing observation to the observation list
        :param observation:
        """
        self.observations.append(observation)

    def shuffle_and_split(self, set_indices: List[int] = None, shuffle: bool = True):
        """
        Returns a new PeptideScores object with a modified observation list
        :param set_indices:
        :param shuffle:
        :return:
        """
        # take the length of the observations
        observations_total: int = len(self.observations)
        # if there are no indices just take the complete set
        if (set_indices is None) or (len(set_indices) == 0):
            set_indices = [0, observations_total]
        # if there is only one entry return a list with only this element
        if len(set_indices) == 1:
            set_indices = [set_indices[0],set_indices[0]+1]
        # when shuffle flag is set to True, prepare an shuffled index array
        shuffled_indices: List[int] = []
        if shuffle:
            shuffled_indices = np.random.permutation(observations_total)
        else:
            shuffled_indices = range(0, observations_total)

        # Puts the elements in a list of lists according to the splits in set_indices
        peptide_score_sets: List[PeptideScores] = []
        for set_size_index in range(0, (len(set_indices) - 1)):
            # calculate the start and end points and remove all boundary issues with modulo
            set_size_start = (set_indices[set_size_index] + observations_total + 1) % (observations_total + 1)
            set_size_end = (set_indices[set_size_index + 1] + observations_total + 1) % (observations_total + 1)
            # create a new PeptideScores object
            peptide_score_sets.append(
                PeptideScores(feature_names=self.feature_names, target_name=self.target_name)
            )
            # fill it with the elements as listed in shuffled_indices
            for i in shuffled_indices[set_size_start:set_size_end]:
                peptide_score_sets[-1].copy_observation(self.observations[i])
        return peptide_score_sets

    def __str__(self):
        return "Peptide Scores for {}:\n{} Observations".format(self.feature_names,len(self.observations))