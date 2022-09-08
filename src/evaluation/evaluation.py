import sys

sys.path.insert(0,"/home/hpc/iwbn/iwbn001h/GenesFRFFinal")


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


if not os.path.exists(config.data_path_to_VIM_matrices + "Weight_Matrix.csv"):
#if 1==1:
    # run GENIE3
    logger.info('Run Genie3')
    start_genie3 = time.time()
    cmd = 'Rscript ' + config.path_to_genie3_R + " " + config.data_path_to_VIM_matrices + " "+ config.path_transcription_factors
    print(cmd)
    os.system(cmd)
    #x = subprocess.check_output(cmd, universal_newlines=True)
    #logger.info('Terminated Genie3 with exit code ', x)
    logger.info('Terminated Genie3 with exit code 1')
    end_genie3 = time.time()

logger.info('Loading Dataset')
data, gene_names, transcription_factors = import_data(config.data_path, config.path_transcription_factors)


if not os.path.exists(os.path.join(config.data_path_to_VIM_matrices, "VIM_federated.csv")):
    # run federated method
    logger.info('Run Federated Approach')
    hospital_data = simulate_different_hospitals(data)
    start_federated = time.time()
    vim_federated = train(hospital_data, gene_names=gene_names, regulators=transcription_factors)
    end_federated = time.time()
    # save VIM federated
    logger.info('saving VIM-matrix from federated approach')
    np.save(os.path.join(config.data_path_to_VIM_matrices, "VIM_federated.npy"), vim_federated)
else:
    # load federated vim matrix
    logger.info('loading VIM matrix from the federated approach')
    path_vim_matrix_federated = os.path.join(config.data_path_to_VIM_matrices, "VIM_federated.npy")
    vim_federated = np.load(path_vim_matrix_federated).astype(float)

logger.info('loading VIM matrix from Genie3')
path = os.path.join(config.data_path_to_VIM_matrices, "Weight_Matrix.csv")
VIM_genie3 = np.loadtxt(path, dtype=str, delimiter=",").astype(float)

logger.info('calculate mse')
mse = mean_squared_error(VIM_genie3, vim_federated)
print("The mse of the two VIM-matrices is: %s" % mse)

logger.info('get linked lists')
edges_genie3 = get_linked_list_federated(VIM_genie3, printing=False)
edges_federated = get_linked_list_federated(vim_federated, printing=False)

logger.info('calculate precision, recall and f1 score')
f1 = [0]
precision = [0]
recall = [0]

num_total = len(edges_federated)
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
        if edges_federated[i] in edges_genie3[:i+1]:
            tp += 1
            fn -= 1
        else:
            fp += 1
        if edges_genie3[i] not in edges_federated[:i+1]:
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

# fig = plt.figure(1)
x = np.arange(0, len(edges_federated)+1)
plt.plot(x, precision, label='precision')
plt.plot(x, recall, label='recall')
plt.plot(x, f1, label='f1')
plt.legend()
plt.xlabel("Number of edges selected")
plt.savefig("precision_recall_f1_scores")
plt.show()

f = open('times.txt', 'w')
f.write("Time Genie3 takes: %s\nTime the federated approach takes: %s" % (
        (end_genie3 - start_genie3), (end_federated - start_federated)))
f.close()
