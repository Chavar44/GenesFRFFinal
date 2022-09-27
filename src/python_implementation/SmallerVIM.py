# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 10:40:52 2022

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
import logging

import src.python_implementation.config as config
import os
from numpy import unravel_index

path = config.data_path_to_VIM_matrices + "VIM_H3.npy"
BigVIM = np.load(path).astype(float)
print("Finished reading Big VIM federated Matrix")
SmallVIM = np.zeros((1637, 60483))

#import regulator names
path_tf=config.path_transcription_factors
tf = np.loadtxt(path_tf, dtype=str)
tf = tf.tolist()
print("Regulators read")

#import TCGA data and read gene_names
path = config.data_path
data = np.loadtxt(path, dtype=str, skiprows=1)
raw_gene_names = data[:, :1]
# Transform raw_gene_names into a list readable by the federated random forest approach
gene_names = []
for xs in raw_gene_names:
    for x in xs:
        # Only reads the first 15 characters in order to correctly compare to the Regulators.txt file
        gene_names.append(x[0:15])
print("Gene names read")



t=0
for ir, row in enumerate(gene_names):
    for j in tf:
        if row == j:
            SmallVIM[t] = BigVIM[ir]
            t = t+1

print("Conversion finished")

path = config.data_path_to_VIM_matrices + "SVIM_H3.npy"

np.save(path, SmallVIM, allow_pickle=False)


