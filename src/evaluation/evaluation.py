
import sys
sys.path.insert(0,"/data_slow/xo53tota/GenesFRFFinal")
from genie3.GENIE3 import *
#from ~/GenesFRFFinal/genie3/GENIE3 import *
from src.python_implementation.main import *
from sklearn.metrics import mean_squared_error
import numpy as np
from matplotlib import pyplot as plt
import time
import logging
import subprocess
import src.python_implementation.config as config
import os

log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)
logger = logging.getLogger(__name__)


logger.info('Loading Dataset')
data, gene_names, transcription_factors = import_data(config.data_path, config.path_transcription_factors)



# run or load federated approach
if not os.path.exists(os.path.join(config.data_path_to_VIM_matrices, "VIM_federated.npy")):
    # run federated method
    logger.info('Run Federated Approach')
    hospital_data = simulate_different_hospitals(data)
    start_federated = time.time()
    vim_federated = train(hospital_data, gene_names=gene_names, regulators=transcription_factors)
    end_federated = time.time()
    # save VIM federated
    logger.info('saving VIM-matrix from federated approach')
    np.save(os.path.join(config.data_path_to_VIM_matrices, "VIM_federated.npy"), vim_federated)
    # save time taken
    f = open(os.path.join(config.path_to_results, 'times.txt'), 'w')
    f.write("Time the federated approach takes: %s\n" % (end_federated - start_federated))
    f.close()
else:
    # load federated vim matrix
    logger.info('loading VIM matrix from the federated approach')
    path_vim_matrix_federated = os.path.join(config.data_path_to_VIM_matrices, "VIM_federated.npy")
    vim_federated = np.load(path_vim_matrix_federated).astype(float)

#Load Genie3Matrix
logger.info('loading VIM matrix from Genie3')
path = os.path.join(config.data_path_to_VIM_matrices, "Weight_Matrix.csv")
VIM_genie3_small = np.loadtxt(path, dtype=str, delimiter=",")
logger.info('loading VIM matrix from Genie3')

#Calculate Genie3 LinkList
VIM_genie3 = np.zeros(vim_federated.shape)
vim_fed_sum = np.sum(vim_federated, axis=1)
sum_reg = 0
j = 0

for i in range(0, len(vim_federated)):
    if vim_fed_sum[i] != 0:
        sum_reg += 1
        VIM_genie3[i] = VIM_genie3_small[j]
        j += 1


file_name_link_list_genie3 = os.path.join(config.data_path_to_VIM_matrices, "LinkListG3.txt")
edges_genie3 = get_linked_list_federated(VIM_genie3, gene_names=gene_names, regulators=transcription_factors,
                                            max_count=config.max_count_link_list,
                                            file_name=file_name_link_list_genie3, printing=False)

del VIM_genie3

# calculate link list from federated approach
logger.info('Calculate linked list from federated approach')
split_name = 'uneven'
if config.split_even:
    split_name = 'even'
file_name_link_list_federated = os.path.join(config.data_path_to_VIM_matrices,
                                             'link_list_' + str(config.number_of_hospitals) + '_' + split_name + ".txt")
edges_federated = get_linked_list_federated(vim_federated, gene_names=gene_names, regulators=transcription_factors,
                                            max_count=config.max_count_link_list,
                                            file_name=file_name_link_list_federated, printing=False)

del data, gene_names, transcription_factors, vim_federated

# TODO auf Genie3 anpassen -> linked list laden (namen der Liste anpassen und richtig laden)
logger.info('loading link list from Genie3')
path = os.path.join(config.data_path_to_VIM_matrices, "Weight_Matrix.csv")
edges_genie3 = np.loadtxt(path, dtype=str, delimiter=",").astype(float)

logger.info('calculate precision, recall and f1 score')
f1 = [0]
precision = [0]
recall = [0]

num_total = min(len(edges_federated), len(edges_genie3))
edges_federated = np.delete(np.asarray(edges_federated), obj=2, axis=1).tolist()
edges_genie3 = np.delete(np.asarray(edges_genie3), obj=2, axis=1).tolist()
tp = 0
tn = 0
fp = 0
fn = 0

for i in range(0, num_total):
    if edges_federated[i] == edges_genie3[i]:
        tp += 1
    else:
        if edges_federated[i] in edges_genie3[:i + 1]:
            tp += 1
            fn -= 1
        else:
            fp += 1
        if edges_genie3[i] not in edges_federated[:i + 1]:
            fn += 1
        else:
            tp += 1
            fp -= 1
    tn = num_total - (tp + fn + fp)
    precision.append(tp / (tp + fp))
    recall.append(tp / (tp + fn))
    if precision[i] + recall[i] != 0:
        f1.append(2 * (precision[i] * recall[i]) / (precision[i] + recall[i]))
    else:
        f1.append(0)

x = np.arange(0, num_total + 1)
plt.plot(x, precision, label='precision')
plt.plot(x, recall, label='recall')
plt.plot(x, f1, label='f1')
plt.legend()
plt.xlabel("Number of edges selected")
file_name_png = 'precision_recall_f1_scores_' + str(config.number_of_hospitals) + '_' + split_name + ".txt"
plt.savefig(os.path.join(config.path_to_results, file_name_png))
plt.show()
