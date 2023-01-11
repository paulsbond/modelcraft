import dataclasses
import gemmi
from ..job import Job


@dataclasses.dataclass
class MapMaskResult:
    mask: gemmi.Ccp4Mask
    seconds: float


class MapMask(Job):
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

    def _result(self) -> MapMaskResult:
        self._check_files_exist("mapmask.mrc")
        return MapMaskResult(
            mask=gemmi.read_ccp4_mask(self._path("mapmask.mrc")),
            seconds=self._seconds,
        )
