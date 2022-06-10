from GENIE3 import *

# data_path = "../data/TCGA-COAD.htseq_fpkm.tsv"
data_path = "data.txt"

# load example data
# data = loadtxt(data_path, dtype=str, skiprows=1)[:, 1:].astype(float)
# f = open(data_path)

data = loadtxt("data.txt" ,skiprows=1)
f = open("data.txt")
gene_names = f.readline()
f.close()
gene_names = gene_names.rstrip('\n').split('\t')[1:]


# run GENIE3
VIM = GENIE3(data)
print(VIM)


# Genes that are used as candidate regulators
# regulators = ['CD19', 'CDH17','RAD51','OSR2','TBX3']
# VIM2 = GENIE3(data,gene_names=gene_names,regulators=regulators)


# Use Extra-Trees method
# tree_method='ET'
# # Number of randomly chosen candidate regulators at each node of a tree
# K = 7
# # Number of trees per ensemble
# ntrees = 50
# # Run the method with these settings
# VIM3 = GENIE3(data,tree_method=tree_method,K=K,ntrees=ntrees)

# obtain more information
help(GENIE3)

# get the predicted ranking of all the regulatory links
get_link_list(VIM)

# show the names of the genes
get_link_list(VIM,gene_names=gene_names)

# # show only the links that are directed from the candidate regulators
# get_link_list(VIM,gene_names=gene_names,regulators=regulators)
#
# # show only the first 5 links only
# get_link_list(VIM,gene_names=gene_names,regulators=regulators,maxcount=5)
#
# # write the predicted links in a file
# get_link_list(VIM,gene_names=gene_names, regulators=regulators,file_name='ranking.txt')

# more info
help(get_link_list)