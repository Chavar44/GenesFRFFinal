from GENIE3 import *
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

logger.info('Run Genie3')
start_federated = time.time()
vim = GENIE3(data, gene_names=gene_names, regulators=transcription_factors, nthreads=4)
end_federated = time.time()

# save VIM federated
logger.info('saving VIM-matrix from federated approach')
np.save(os.path.join(config.data_path_to_VIM_matrices, "VIM_Genie3.npy"), vim, allow_pickle=False)

# datetime object containing current date and time
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
f = open('times.txt', 'a')
f.write("%s Time Genie3 takes: %s" % (dt_string, (end_federated - start_federated)))
f.close()