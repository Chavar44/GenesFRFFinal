# imports

# data_path = "../../data/TCGA-COAD.htseq_fpkm.tsv"
# data_path = "/genie3/data.txt"
data_path = "/data_slow/xo53tota/GenesFRFFinal/data/NewDataset.tsv"
path_transcription_factors = '/data_slow/xo53tota/GenesFRFFinal/data/Regulators.txt'
path_to_genie3_R = '/data_slow/xo53tota/GenesFRFFinal/genie3/Genie3.R'
data_path_to_VIM_matrices = "/data_slow/xo53tota/GenesFRFFinal/data/"
path_to_results = '/data_slow/xo53tota/GenesFRFFinal/results/'

number_trees = 500
tree_method = "RF"

number_of_hospitals = 10
split_even = True
# is the list of how to split uneven data: len must be of number_of_hospitals
split_uneven = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
#For reminder purposes split_uneven = [0.1, 0.075, 0.075, 0.225, 0.125, 0.075, 0.075, 0.1, 0.075, 0.075]
#For reminder purposes split_uneven = [0.1, 0.2, 0.7]


max_count_link_list = 1000000

#19574 coding genes in dataset
#density = Edges_analysed/(Genes*regulators)
density = (max_count_link_list/(19574*1637)) * 100