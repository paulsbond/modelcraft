import dataclasses
import os
import gemmi

from .nucleofind import NucleoFindResult
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


@dataclasses.dataclass
class CootNucleoFindRSRResult:
    structure: gemmi.Structure
    seconds: int


class CootNucleoFindRSR(Job):
    def __init__(self,
                 fsigf: DataItem,
                 phases: DataItem,
                 fphi: DataItem = None,
                 freer: DataItem = None,
                 structure: gemmi.Structure = None,
                 nucleofind_result: NucleoFindResult = None,
                 chain: str = "",
                 res_range_start: int = None,
                 res_range_end: int = None):
        super().__init__("coot-mini-rsr")

        self.fsigf = fsigf
        self.phases = phases
        self.fphi = fphi
        self.freer = freer
        self.structure = structure
        self.chain = chain
        self.res_range_start = res_range_start
        self.res_range_end = res_range_end
        self.nucleofind_result = nucleofind_result
        self.other_chains = []

    def _setup(self):
        data_items = [self.fsigf, self.phases, self.fphi, self.freer]
        write_mtz(self._path("hklin.mtz"), data_items)

        if self.nucleofind_result:
            mtz = gemmi.read_mtz_file(self._path("hklin.mtz"))
            grid = mtz.transform_f_phi_to_map(self.fphi.label(0), self.fphi.label(1))
            grid.normalize()
            grid = self._weight_map(grid)
            grid_name = "grid.map"
            self._save_map(grid, grid_name)
            self._args += ["--mapin", grid_name]
        else:
            if self.fphi is not None:
                self._args += ["--f", self.fphi.label(0)]
                self._args += ["--phi", self.fphi.label(1)]
            self._args += ["--hklin", "hklin.mtz"]

        if self.structure is not None:
            write_mmcif(self._path("xyzin.cif"), self.structure)
            self._args += ["--pdbin", "xyzin.cif"]

        self._args += ["--pdbout", "xyzout.cif"]
        if self.chain:
            self._args += ["--chain-id", self.chain]
        if self.res_range_start is not None:
            self._args += ["--resno-start", self.res_range_start]
        if self.res_range_end is not None:
            self._args += ["--resno-end", self.res_range_end]

    def _weight_map(self, grid: gemmi.FloatGrid) -> gemmi.FloatGrid:
        for point in grid:
            position = grid.point_to_position(point)
            phosphate_value = self.nucleofind_result.predicted_phosphate_map.grid.interpolate_value(position)
            sugar_value = self.nucleofind_result.predicted_sugar_map.grid.interpolate_value(position)
            base_value = self.nucleofind_result.predicted_base_map.grid.interpolate_value(position)
            point.value = point.value * (sugar_value + phosphate_value + base_value)
        return grid

    def _save_map(self, grid: gemmi.FloatGrid, name: str = "grid.map"):
        m = gemmi.Ccp4Map()
        m.grid = grid
        m.update_ccp4_header()
        m.write_ccp4_map(self._path(name))


    def _result(self) -> gemmi.Structure:
        self._check_files_exist("xyzout.cif")

        refined_structure = read_structure(self._path("xyzout.cif"))

        refined_chain = refined_structure[0].find_chain(self.chain)
        original_chain = self.structure[0].find_chain(self.chain)
        superposition = gemmi.calculate_superposition(original_chain.whole(), refined_chain.whole(),
                                                      gemmi.PolymerType.Rna, gemmi.SupSelect.All, trim_cycles=5)
        threshold = 2
        if superposition.rmsd > threshold:
            return CootNucleoFindRSR(self.structure, self._seconds)
        return CootNucleoFindRSR(refined_structure, self._seconds)
