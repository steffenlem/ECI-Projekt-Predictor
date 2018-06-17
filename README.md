# ECI-Projekt-Predictor
Epitope Binding Predictor

This epitope binding predictor uses an ANN to predict whether a peptide binds to a
MHC-I or not.

arguments
--train <filepath>  Path to the training file
--input <filepath>  Path to the input file
--output <filepath> Path to the output file

Example for running with the test_input_easy.txt:
Go to "Run > Run...", "Edit Configurations...", then add this to the "Parameters" field:
--input "data\test_input_easy.txt"