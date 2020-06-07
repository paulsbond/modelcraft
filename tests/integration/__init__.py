import os


def ccp4_path(*paths: str) -> str:
    if "CCP4" not in os.environ:
        raise EnvironmentError("CCP4 environment not set")
    return os.path.join(os.environ["CCP4"], *paths)
