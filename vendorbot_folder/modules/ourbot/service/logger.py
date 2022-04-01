import sys
import logging

logging.basicConfig(
    # filename='my_runtime_log.log', # saving log to filename
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO  # DEBUG
)
logging.info('logger started')
log = logging.getLogger()
