from genie3.GENIE3 import *
from src.python_implementation.main import *
import numpy as np

data_path = "/media/sf_Projekt_BIONETS/federated-inference-of-grns/genie3/data.txt"
data = import_data(data_path)

# run GENIE3
VIM_genie3 = GENIE3(data)

number_patients = data.shape[0]
number_genes = data.shape[1]
hospital_data = simulate_different_hospitals(data)
vim_federated = train(hospital_data, number_genes)

print(VIM_genie3 == vim_federated)
print(np.abs(VIM_genie3-vim_federated))

get_link_list(VIM_genie3)
print("federated")
get_linked_list_federated(vim_federated)
