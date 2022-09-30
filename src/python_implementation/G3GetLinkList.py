# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 10:38:09 2022

@author: chvr4
"""
import sys
sys.path.insert(0,"/data_slow/xo53tota/GenesFRFFinal")
from genie3.GENIE3 import *
#from ~/GenesFRFFinal/genie3/GENIE3 import *
from src.python_implementation.main import *

import numpy as np

import time
import logging

import src.python_implementation.config as config

from numpy import unravel_index

max_count_link_list = 5000000



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

#Get order of regulators in matrix
Regulators = []
j=0
for i in gene_names:
    for j in tf:
        if i == j:
            Regulators.append(j)

print("Regulators in order") 


#import G3 VIM matrix
path = config.data_path_to_VIM_matrices + "VIM_Federated_10_even.npy"
G3 = np.load(path).astype(float)
print("Finished reading F Matrix")

path = config.data_path_to_VIM_matrices + "Federated_10_even_linked_list_5000000.txt"
with open(path, 'w') as f:
    
    for i in range(0,max_count_link_list):
        coords = unravel_index(G3.argmax(), G3.shape) #get coordinates for the VIM Matrix in order to know the exact gene
        line =  Regulators[coords[0]]  + " " + gene_names[coords[1]] + " " + str(G3[coords[0]][coords[1]])
        f.write(line)
        f.write('\n')

        print(str(i) + " " + line)
        print(coords[0])
        print(coords[0])
        G3[coords[0]][coords[1]] = 0
    
f.close()



