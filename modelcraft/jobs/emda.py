import dataclasses
import gemmi
from ..job import Job
from ..maps import read_map


@dataclasses.dataclass
class EmdaMapMaskResult:
    mask: gemmi.Ccp4Map
    seconds: float


class EmdaMapMask(Job):
    def __init__(
        self,
        density: gemmi.Ccp4Map,
        kernel_radius: int = 10,
        resolution: float = 10.0,
    ):
        super().__init__("emda")
        self.density = density
        self.kernel_radius = kernel_radius
        self.resolution = resolution

    def _setup(self) -> None:
        self.density.write_ccp4_map(self._path("map.cpp4"))
        self._args += ["mapmask"]
        self._args += ["--map", "map.cpp4"]
        self._args += ["--knl", str(self.kernel_radius)]
        self._args += ["--res", str(self.resolution)]

    def _result(self) -> EmdaMapMaskResult:
        self._check_files_exist("mapmask.mrc")
        return EmdaMapMaskResult(
            mask=read_map(self._path("mapmask.mrc")),
            seconds=self._seconds,
        )
