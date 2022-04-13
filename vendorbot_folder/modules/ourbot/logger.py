import logging
import sys
log_level = logging.INFO


class Logger(object):
    def __init__(self):
        log = logging.getLogger('')
        log.setLevel(log_level)

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        console = sys.stderr
        if console is not None:
            console = logging.StreamHandler(console)
            console.setLevel(log_level)
            console.setFormatter(formatter)
            log.addHandler(console)

Logger()

log = logger = logging.getLogger('')
