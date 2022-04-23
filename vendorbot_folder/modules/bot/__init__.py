import os


def run_async():
    return os.getenv("RUN_ASYNC", "true") == "true"
