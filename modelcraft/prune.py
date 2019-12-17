from modelcraft.job import Job
import os
import shutil


class Prune(Job):
    def __init__(self, directory, xyzin, hklin, chains_only=False):
        super().__init__(directory)
        self.xyzout = self.path("xyzout.pdb")
        path = os.path.join(os.path.dirname(__file__), "coot_prune.py")
        with open(path) as coot_prune:
            script = coot_prune.read()
        if chains_only:
            script += "deleted = prune(0, 1, 2, residues=False, sidechains=False)\n"
        else:
            script += "deleted = prune(0, 1, 2)\n"
        script += "write_pdb_file(0, '%s')\n" % self.xyzout
        script += "exit()"
        with open(self.path("script.py"), "w") as f:
            f.write(script)
        arguments = [
            xyzin,
            hklin.path,
            "--no-graphics",
            "--no-guano",
            "--no-state-script",
            "--script", self.path("script.py")
        ]
        self.run("coot", arguments)
        shutil.rmtree("coot-backup", ignore_errors=True)
        shutil.rmtree("coot-download", ignore_errors=True)
