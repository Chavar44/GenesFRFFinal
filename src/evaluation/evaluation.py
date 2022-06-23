from genie3.GENIE3 import *
from src.python_implementation.main import *
from sklearn.metrics import *
import numpy as np

data_path = "/media/sf_Projekt_BIONETS/federated-inference-of-grns/genie3/data.txt"
data = import_data(data_path)
threshold = 0.01

# run GENIE3
VIM_genie3 = GENIE3(data)

# run federated method
number_patients = data.shape[0]
number_genes = data.shape[1]
hospital_data = simulate_different_hospitals(data)
vim_federated = train(hospital_data, number_genes,number_patients)

# get diff of VIM's
diff_of_VIM = VIM_genie3 - vim_federated
# if diff is greater than threshold, it is a "false" prediction
# if 0 the prediction is true, else 1
diff_of_VIM = np.where(diff_of_VIM < threshold, 0, 1)
print(diff_of_VIM)
print("Number of false predicted %s, in percentage: %s" % (np.sum(diff_of_VIM), np.sum(diff_of_VIM)/(shape(diff_of_VIM)[0]**2 - shape(diff_of_VIM)[0])*100))
mse = mean_squared_error(VIM_genie3, vim_federated)
print(mse)

# get_link_list(VIM_genie3)
# print("federated")
# get_linked_list_federated(vim_federated)
