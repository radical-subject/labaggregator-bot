import os
import sys
from logging import getLogger, StreamHandler, Formatter, INFO, DEBUG

log_level = DEBUG


class Logger(object):
    def __init__(self):
        log = getLogger('')

        formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # Use console for log output
        if os.getenv('NOCONSOLE') == '1':
            console = sys.stderr
            if console is not None:
                # Logging to console and file both
                console = StreamHandler(console)
                console.setLevel(log_level)
                console.setFormatter(formatter)
                log.addHandler(console)


log = logger = getLogger('')
