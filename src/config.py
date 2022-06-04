# imports

data_path = "../data/TCGA-COAD.htseq_fpkm.tsv"

number_trees = 500
tree_method = "RF"

number_of_hospitals = 3
split_even = True
# is the list of how to split uneven data: len must be of number_of_hospitals
split_uneven = [0.1, 0.2, 0.7]
