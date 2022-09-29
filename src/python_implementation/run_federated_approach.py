import sys
sys.path.insert(0,"/data_slow/xo53tota/GenesFRFFinal")
from src.python_implementation.main import *
import numpy as np
import time
import logging
import src.python_implementation.config as config
import os
from datetime import datetime

log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)
logger = logging.getLogger(__name__)


logger.info('Loading Dataset')
data, gene_names, transcription_factors = import_data(config.data_path, config.path_transcription_factors)

# run federated method
logger.info('Run Federated Approach')
hospital_data = simulate_different_hospitals(data)
start_federated = time.time()
vim_federated = train(hospital_data, gene_names=gene_names, regulators=transcription_factors, parallelize_hospitals=config.number_of_hospitals)
end_federated = time.time()

# save VIM federated
logger.info('saving VIM-matrix from federated approach')
np.save(os.path.join(config.data_path_to_VIM_matrices, "VIM_federated.npy"), vim_federated, allow_pickle=False)

# datetime object containing current date and time
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
f = open('times.txt', 'a')
f.write("%s Time Federated Approach takes: %s" % (dt_string, (end_federated - start_federated)))
f.close()