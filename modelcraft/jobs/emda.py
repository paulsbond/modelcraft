import dataclasses
import gemmi
from ..job import Job


@dataclasses.dataclass
class MapMaskResult:
    mask: gemmi.FloatGrid
    seconds: float


class MapMask(Job):
    def __init__(
        self,
        density: gemmi.FloatGrid,
        kernel_radius: int = 10,
        resolution: float = 10.0,
    ):
        super().__init__("emda")
        self.density = density
        self.kernel_radius = kernel_radius
        self.resolution = resolution

    def _setup(self) -> None:
        _write_map(self.density, self._path("mapin.cpp4"))
        self._args += ["mapmask"]
        self._args += ["--map", "mapin.cpp4"]
        self._args += ["--knl", str(self.kernel_radius)]
        self._args += ["--res", str(self.resolution)]

    def _result(self) -> MapMaskResult:
        self._check_files_exist("mapmask.mrc")
        return MapMaskResult(
            mask=gemmi.read_ccp4_map(self._path("mapmask.mrc")).grid,
            seconds=self._seconds,
        )


def _write_map(grid: gemmi.FloatGrid, path: str):
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(path)
