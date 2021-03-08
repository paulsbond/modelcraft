import os


def setup_environ():
    for variable in ("CCP4", "CLIBD"):
        if variable not in os.environ:
            raise EnvironmentError(variable + " environment variable not set")
    os.environ["COOT_N_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
