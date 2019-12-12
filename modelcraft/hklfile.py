import os


class HklFile():
    def __init__(self, path, fsigf, abcd=None, phifom=None, fphi=None):
        self.path = os.path.abspath(path)
        self.fsigf = fsigf
        self.abcd = abcd
        self.phifom = phifom
        self.fphi = fphi
