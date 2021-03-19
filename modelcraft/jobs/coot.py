import dataclasses
import os
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class CootResult:
    structure: gemmi.Structure


class Coot(Job):
    def __init__(
        self,
        structure: gemmi.Structure,
        fphi_best: DataItem,
        fphi_diff: DataItem,
        script: str,
    ):
        super().__init__("coot")
        self.structure = structure
        self.fphi_best = fphi_best
        self.fphi_diff = fphi_diff
        self.script = script

    def _setup(self) -> None:
        script_lines = [
            "handle_read_draw_molecule('xyzin.cif')\n",
            "make_and_draw_map('hklin.mtz',",
            f"    '{self.fphi_best.label(0)}', '{self.fphi_best.label(1)}', '', 0, 0)\n",
            "make_and_draw_map('hklin.mtz',",
            f"    '{self.fphi_diff.label(0)}', '{self.fphi_diff.label(1)}', '', 0, 1)\n",
            "turn_off_backup(0)\n",
            "try:\n",
        ]
        for line in self.script.split("\n"):
            script_lines.append("    %s\n" % line)
        script_lines += [
            "    write_cif_file(0, 'xyzout.cif')\n",
            "except:\n",
            "    import traceback\n",
            "    traceback.print_exc()\n",
            "    coot_real_exit(1)\n",
            "coot_real_exit(0)\n",
        ]
        with open(self._path("script.py"), "w") as script_file:
            script_file.writelines(script_lines)
        write_mtz(self._path("hklin.mtz"), [self.fphi_best, self.fphi_diff])
        write_mmcif(self._path("xyzin.cif"), self.structure)
        self._args += ["--no-graphics"]
        self._args += ["--no-guano"]
        self._args += ["--no-state-script"]
        self._args += ["--script", "script.py"]

    def _result(self) -> CootResult:
        return CootResult(structure=read_structure(self._path("xyzout.cif")))


class Prune(Coot):
    def __init__(
        self,
        structure: gemmi.Structure,
        fphi_best: DataItem,
        fphi_diff: DataItem,
        chains_only: bool = False,
    ):
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "prune.py")
        with open(path) as stream:
            script = stream.read()
        if chains_only:
            script += "prune(0, 1, 2, residues=False, sidechains=False)\n"
        else:
            script += "prune(0, 1, 2)\n"
        super().__init__(structure, fphi_best, fphi_diff, script)


class FixSideChains(Coot):
    def __init__(
        self,
        structure: gemmi.Structure,
        fphi_best: DataItem,
        fphi_diff: DataItem,
    ):
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "prune.py")
        with open(path) as stream:
            script = stream.read()
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "sidechains.py")
        with open(path) as stream:
            script += "\n\n%s\n" % stream.read()
        script += "fix_side_chains(0, 1, 2)\n"
        super().__init__(structure, fphi_best, fphi_diff, script)
