import dataclasses
import os
import gemmi
from ..job import Job
from ..reflections import DataItem, write_mtz
from ..structure import read_structure, write_mmcif


@dataclasses.dataclass
class CootResult:
    structure: gemmi.Structure
    seconds: float


class Coot(Job):
    def __init__(self, script: str, structures: list, fphis: list):
        super().__init__("coot")
        self.script = script
        self.structures = structures
        self.fphis = fphis

    def _setup(self) -> None:
        script_lines = [
            "try:\n",
            "    COOT1 = True\n",
            "    try:\n",
            "        import coot_utils\n",
            "    except NameError:\n",
            "        COOT1 = False\n",
            "    if COOT1:\n",
            "        from coot import *\n",
            "    turn_off_backup(0)\n",
        ]
        for i, structure in enumerate(self.structures):
            write_mmcif(self._path(f"xyzin{i}.cif"), structure)
            script_lines += [
                f"    IMOL{i} = handle_read_draw_molecule('xyzin{i}.cif')\n"
            ]
        for i, fphi in enumerate(self.fphis):
            write_mtz(self._path(f"hklin{i}.mtz"), [fphi])
            script_lines += [
                f"    IMAP{i} = make_and_draw_map('hklin{i}.mtz', "
                f"'{fphi.label(0)}', '{fphi.label(1)}', '', 0, 0)\n"
            ]
        script_lines += ["    set_imol_refinement_map(IMAP0)\n"]
        for line in self.script.split("\n"):
            script_lines += [f"    {line}\n"]
        script_lines += [
            "    write_cif_file(IMOL0, 'xyzout.cif')\n",
            "    coot_real_exit(0)\n",
            "except:\n",
            "    import traceback\n",
            "    traceback.print_exc()\n",
            "    coot_real_exit(1)\n",
        ]
        with open(self._path("script.py"), "w") as script_file:
            script_file.writelines(script_lines)
        self._args += ["--no-graphics"]
        self._args += ["--no-guano"]
        self._args += ["--no-state-script"]
        self._args += ["--script", "script.py"]

    def _result(self) -> CootResult:
        self._check_files_exist("xyzout.cif")
        return CootResult(
            structure=read_structure(self._path("xyzout.cif")),
            seconds=self._seconds,
        )


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
            script += "prune(IMOL0, IMAP0, IMAP1, residues=False, sidechains=False)\n"
        else:
            script += "prune(IMOL0, IMAP0, IMAP1)\n"
        super().__init__(
            script=script, structures=[structure], fphis=[fphi_best, fphi_diff]
        )


class FixSideChains(Coot):
    def __init__(
        self, structure: gemmi.Structure, fphi_best: DataItem, fphi_diff: DataItem
    ):
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "prune.py")
        with open(path) as stream:
            script = stream.read()
        path = os.path.join(os.path.dirname(__file__), "..", "coot", "sidechains.py")
        with open(path) as stream:
            script += "\n\n%s\n" % stream.read()
        script += "fix_side_chains(IMOL0, IMAP0, IMAP1)\n"
        super().__init__(
            script=script, structures=[structure], fphis=[fphi_best, fphi_diff]
        )
