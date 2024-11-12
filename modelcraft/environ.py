import os


def setup_environ():
    for variable in ("CCP4", "CLIB", "CLIBD_MON"):
        if variable not in os.environ:
            raise EnvironmentError(variable + " environment variable not set")
    os.environ["LD_LIBRARY_PATH"] = os.environ["CLIB"]
