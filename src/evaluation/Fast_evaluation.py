# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 18:49:56 2022

@author: chvr4
"""

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

split_name = 'uneven'
if config.split_even:
    split_name = 'even'

logger.info('loading link list from Federated Approach')
path = os.path.join(config.data_path_to_VIM_matrices, "F_uneven_3_linked_list_1100000.txt")
edges_federated = np.loadtxt(path, dtype=str, delimiter=" ")

logger.info('loading link list from Genie3')
path = os.path.join(config.data_path_to_VIM_matrices, "G3_linked_list_1500000.txt")
edges_genie3 = np.loadtxt(path, dtype=str, delimiter=" ")


logger.info('calculate precision, recall and f1 score')
f1 = [0]
precision = [0]
recall = [0]

num_total = config.max_count_link_list

edges_federated = np.delete(np.asarray(edges_federated), obj=2, axis=1).tolist()
edges_genie3 = np.delete(np.asarray(edges_genie3), obj=2, axis=1).tolist()


tp = 0
tn = 0
fp = 0
fn = 0

for i in range(0, num_total):
    print(i)
    #print("Federated: " + edges_federated[i] + "Genie3: " +  edges_genie3[i] )
    if edges_federated[i] == edges_genie3[i]:
        tp += 1
    else:
        if edges_federated[i] in edges_genie3[:i + 1]:
            tp += 1
        else:
            fp += 1
        if edges_genie3[i] in edges_federated[:i + 1]:
            tp += 1
            fn -= 1
            fp -= 1
        else:
            fn += 1
    tn = num_total - (tp + fn + fp)
    print("tp = " + str(tp))
    print("fp = " + str(fp))
    print("fn = " + str(fn))

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
plt.grid()
plt.xlabel("Number of edges selected")
plt.title("Density(%) = " + "%.4f" % config.density)
file_name_png = 'eval_' + str(config.number_of_hospitals) + '_' + split_name + '_' + str(config.max_count_link_list) + ".png"
plt.savefig(os.path.join(config.path_to_results, file_name_png))
plt.show()
