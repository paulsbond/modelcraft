import os


def setup_environ():
    for variable in ("CCP4", "CLIB", "CLIBD_MON"):
        if variable not in os.environ:
            raise EnvironmentError(variable + " environment variable not set")
    os.environ["LD_LIBRARY_PATH"] = os.environ["CLIB"]
    os.environ["COOT_N_THREADS"] = "1"
    os.environ["GOTO_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
