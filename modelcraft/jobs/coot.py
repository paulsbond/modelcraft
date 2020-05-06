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
        xyzin = self.path("xyzin.pdb")  # TODO: Change to CIF
        hklin = self.path("hklin.mtz")
        script_path = self.path("script.py")
        xyzout = self.path("xyzout.pdb")  # TODO: Change to CIF
        script = (
            f"read_pdb('{xyzin}')\n"
            f"make_and_draw_map('{hklin}', '{fphi_best.label(0)}', '{fphi_best.label(1)}', '', 0, 0)\n"
            f"make_and_draw_map('{hklin}', '{fphi_diff.label(0)}', '{fphi_diff.label(1)}', '', 0, 1)\n"
            f"{script}\n"
            f"write_pdb_file(0, '{xyzout}')\n"  # TODO: write_cif_file
            "exit()\n"
        )

        with open(script_path, "w") as script_file:
            script_file.write(script)
        structure.write_pdb(xyzin)  # TODO: Change to CIF
        # write_mmcif(xyzin, structure)
        write_mtz(hklin, [fphi_best, fphi_diff])

        args = []
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


class FixSideChains(Coot):
    def __init__(self, structure: gemmi.Structure, fphi_best: FPhi, fphi_diff: FPhi):
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "prune.py")
        with open(path) as script_file:
            script = script_file.read()
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "sidechains.py")
        with open(path) as script_file:
            script += "\n\n%s\n" % script_file.read()
        script += "fix_side_chains(0, 1, 2)\n"
        super().__init__(structure, fphi_best, fphi_diff, script)
