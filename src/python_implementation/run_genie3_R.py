import time
import logging
import subprocess
import src.python_implementation.config as config
from datetime import datetime

log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=log_fmt)
logger = logging.getLogger(__name__)


logger.info('Loading Dataset')

# run GENIE3
logger.info('Run Genie3')
start_genie3 = time.time()
cmd = ['Rscript', config.path_to_genie3_R, config.data_path, config.path_transcription_factors]
x = subprocess.check_output(cmd, universal_newlines=True)
logger.info('Terminated Genie3 with exit code %s', x)
end_genie3 = time.time()

# datetime object containing current date and time
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
f = open('times.txt', 'a')
f.write("%s Time Genie3 takes: %s" % (dt_string, (end_genie3 - start_genie3)))
f.close()
