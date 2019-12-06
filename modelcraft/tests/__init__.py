import os


def data_path(filename):
    return os.path.join(os.path.dirname(__file__), "data", filename)
