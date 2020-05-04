import os
import shutil
import gemmi
from ..structure import write_mmcif
from ..reflections import FPhi, write_mtz
from .job import Job


class Coot(Job):
    def __init__(
        self, structure: gemmi.Structure, fphi_best: FPhi, fphi_diff: FPhi, script: str
    ):
        super().__init__()
        xyzin = self.path("xyzin.cif")
        hklin = self.path("hklin.mtz")
        script_path = self.path("script.py")
        xyzout = self.path("xyzout.cif")

        script += f"\nwrite_pdb_file(0, '{xyzout}')\n"
        script += "exit()"
        with open(script_path, "w") as script_file:
            script_file.write(script)
        write_mmcif(xyzin, structure)
        write_mtz(hklin, [fphi_best, fphi_diff])

        args = []
        # xyzin is not being read (maybe due to cif)
        # change both this and the map reading to be done with scripting commands)
        args += [xyzin]
        args += [hklin]
        args += ["--no-graphics"]
        args += ["--no-guano"]
        args += ["--no-state-script"]
        args += ["--script", script_path]
        self.run("coot", args)
        self.structure = gemmi.read_structure(xyzout)
        shutil.rmtree("coot-backup", ignore_errors=True)
        shutil.rmtree("coot-download", ignore_errors=True)
        self.finish()


class Prune(Coot):
    def __init__(
        self,
        structure: gemmi.Structure,
        fphi_best: FPhi,
        fphi_diff: FPhi,
        chains_only: bool = False,
    ):
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "prune.py")
        with open(path) as script_file:
            script = script_file.read()
        if chains_only:
            script += "prune(0, 1, 2, residues=False, sidechains=False)\n"
        else:
            script += "prune(0, 1, 2)\n"
        super().__init__(structure, fphi_best, fphi_diff, script)
