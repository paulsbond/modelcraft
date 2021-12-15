import dataclasses
from ..job import Job
from ..reflections import DataItem, write_mtz


@dataclasses.dataclass
class PhaseMatchResult:
    f_map_correlation: float
    seconds: float


class PhaseMatch(Job):
    def __init__(self, fsigf: DataItem, phases1: DataItem, phases2: DataItem):
        super().__init__("cphasematch")
        self.fsigf = fsigf
        self.phases1 = phases1
        self.phases2 = phases2

    def _setup(self) -> None:
        self._args += ["-mtzin", "hklin.mtz"]
        self._args += ["-colin-fo", "F,SIGF"]
        arg1 = "-colin-hl-1" if self.phases1.types == "AAAA" else "-colin-phifom-1"
        arg2 = "-colin-hl-2" if self.phases2.types == "AAAA" else "-colin-phifom-2"
        label1 = "HLA1,HLB1,HLC1,HLD1" if self.phases1.types == "AAAA" else "PHI1,FOM1"
        label2 = "HLA2,HLB2,HLC2,HLD2" if self.phases2.types == "AAAA" else "PHI2,FOM2"
        self._args += [arg1, label1]
        self._args += [arg2, label2]
        self._args += ["-mtzout", "hklout.mtz"]
        write_mtz(
            path=self._path("hklin.mtz"),
            items=[self.fsigf, self.phases1, self.phases2],
            labels=["F,SIGF", label1, label2],
        )

    def _result(self) -> PhaseMatchResult:
        self._check_files_exist("hklout.mtz")
        path = self._path("stdout.txt")
        with open(path) as stream:
            for line in stream:
                if line[:19] == "Overall statistics:":
                    keys = next(stream).split()
                    values = next(stream).split()
                    stats = dict(zip(keys, values))
                    return PhaseMatchResult(
                        f_map_correlation=float(stats["wFcorr"]),
                        seconds=self._seconds,
                    )
        raise RuntimeError(f"Could not find overall statistics in {path}")
