# imports

# data_path = "../../data/TCGA-COAD.htseq_fpkm.tsv"
# data_path = "/genie3/data.txt"
data_path = "/media/sf_Projekt_BIONETS/federated-inference-of-grns/data/TCGA-COAD.htseq_fpkm.tsv"
path_transcription_factors = '/media/sf_Projekt_BIONETS/federated-inference-of-grns/data/Regulators.txt'
path_to_genie3_R = '/media/sf_Projekt_BIONETS/federated-inference-of-grns/genie3/Genie3.R'
data_path_to_VIM_matrices = "/media/sf_Projekt_BIONETS/federated-inference-of-grns/data/"

number_trees = 500
tree_method = "RF"

number_of_hospitals = 3
split_even = True
# is the list of how to split uneven data: len must be of number_of_hospitals
split_uneven = [0.1, 0.2, 0.7]
