import os
import gemmi
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif
from .job import Job


class Coot(Job):
    def __init__(
        self,
        name: str,
        structure: gemmi.Structure,
        fphi_best: DataItem,
        fphi_diff: DataItem,
        script: str,
    ):
        super().__init__(name)
        xyzin = self.path("xyzin.cif")
        hklin = self.path("hklin.mtz")
        script_path = self.path("script.py")
        xyzout = self.path("xyzout.cif")
        script_lines = [
            f"handle_read_draw_molecule('{xyzin}')\n",
            f"make_and_draw_map('{hklin}',",
            f"    '{fphi_best.label(0)}', '{fphi_best.label(1)}', '', 0, 0)\n",
            f"make_and_draw_map('{hklin}',",
            f"    '{fphi_diff.label(0)}', '{fphi_diff.label(1)}', '', 0, 1)\n",
            "turn_off_backup(0)\n",
            "try:\n",
        ]
        for line in script.split("\n"):
            script_lines.append("    %s\n" % line)
        script_lines += [
            f"    write_cif_file(0, '{xyzout}')\n",
            "except:\n",
            "    import traceback\n",
            "    traceback.print_exc()\n",
            "    coot_real_exit(1)\n",
            "coot_real_exit(0)\n",
        ]

        with open(script_path, "w") as script_file:
            script_file.writelines(script_lines)
        write_mmcif(xyzin, structure)
        write_mtz(hklin, [fphi_best, fphi_diff])

        args = []
        args += ["--no-graphics"]
        args += ["--no-guano"]
        args += ["--no-state-script"]
        args += ["--script", script_path]
        self.run("coot", args)
        self.structure = read_structure(xyzout)
        self.finish()


class Prune(Coot):
    def __init__(
        self,
        structure: gemmi.Structure,
        fphi_best: DataItem,
        fphi_diff: DataItem,
        chains_only: bool = False,
    ):
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "prune.py")
        with open(path) as script_file:
            script = script_file.read()
        if chains_only:
            name = "prune_chains"
            script += "prune(0, 1, 2, residues=False, sidechains=False)\n"
        else:
            name = "prune"
            script += "prune(0, 1, 2)\n"
        super().__init__(name, structure, fphi_best, fphi_diff, script)


class FixSideChains(Coot):
    def __init__(
        self, structure: gemmi.Structure, fphi_best: DataItem, fphi_diff: DataItem
    ):
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "prune.py")
        with open(path) as script_file:
            script = script_file.read()
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "sidechains.py")
        with open(path) as script_file:
            script += "\n\n%s\n" % script_file.read()
        script += "fix_side_chains(0, 1, 2)\n"
        super().__init__("fix_side_chains", structure, fphi_best, fphi_diff, script)
