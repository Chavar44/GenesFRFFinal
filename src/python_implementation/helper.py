import numpy as np
import os
import src.python_implementation.config as config

path_vim_matrix_federated = os.path.join(config.data_path_to_VIM_matrices, "VIM_federated.npy")
vim_federated = np.load(path_vim_matrix_federated).astype(float)

path = os.path.join(config.data_path_to_VIM_matrices, "Weight_Matrix.csv")
VIM_genie3_small = np.loadtxt(path, dtype=str, delimiter=",").astype(float)

VIM_genie3 = np.zeros(vim_federated.shape)
vim_fed_sum = np.sum(vim_federated, axis=1)
sum_reg = 0
j = 0

for i in range(0, len(vim_federated)):
    if vim_fed_sum[i] != 0:
        sum_reg += 1
        VIM_genie3[i] = VIM_genie3_small[j]
        j += 1

np.save(os.path.join(config.data_path_to_VIM_matrices, "VIM_Genie3.npy"), vim_federated, allow_pickle=False)
