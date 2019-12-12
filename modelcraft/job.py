import os


class Job:
    def __init__(self, directory):
        os.makedirs(directory, exist_ok=True)
        self.directory = os.path.abspath(directory)
        self.stdout = self.path("stdout.txt")
        self.stderr = self.path("stderr.txt")
        self.xyzout = self.path("xyzout.pdb")

    def path(self, *paths):
        return os.path.join(self.directory, *paths)
