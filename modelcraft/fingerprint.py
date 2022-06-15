from math import cos, floor, radians, sin
from typing import Iterator, List
import gemmi
from numpy import arange
from scipy.spatial.transform import Rotation


class Fingerprint:
    def __init__(self, high: List[gemmi.Position], low: List[gemmi.Position]):
        self._high = high
        self._low = low
        self._grid = None
        self._grid_scores = []

    def set_grid(self, grid: gemmi.FloatGrid) -> None:
        self._grid = grid
        self._grid_scores = []
        for x in range(0, 1, 1 / 45):
            for y in range(0, 1, 1 / 45):
                for z in range(0, 1, 1 / 45):
                    position = grid.unit_cell.orthogonalize(gemmi.Fractional(x, y, z))
                    transform = gemmi.Transform(gemmi.Mat33(), position)
                    self._grid_scores.append(self.score(transform))
        self._grid_scores.sort()

    def radius(self) -> float:
        return max(pos.length() for pos in self._high + self._low) + 1.0

    def high_values(self, transform: gemmi.Transform) -> Iterator[float]:
        for position in self._high:
            transformed = gemmi.Position(*transform.apply(position))
            yield self._grid.tricubic_interpolation(transformed)

    def low_values(self, transform: gemmi.Transform) -> Iterator[float]:
        for position in self._low:
            transformed = gemmi.Position(*transform.apply(position))
            yield self._grid.tricubic_interpolation(transformed)

    def score(self, transform: gemmi.Transform) -> float:
        return sum(self.high_values(transform)) - sum(self.low_values(transform))

    def search(self) -> List[gemmi.Transform]:
        transforms = []
        value_cutoff = self._grid.array.mean() + self._grid.array.std()
        for point in self._grid.masked_asu():
            if point.value < value_cutoff:
                continue
            position = self._grid.point_to_position(point)
            best_score = 0
            best_transform = None
            for rotation in _rotations(step=18):
                transform = gemmi.Transform(rotation, position)
                high_values = self.high_values(transform)
                low_values = self.low_values(transform)
                min_high = max_low = 0
                for high_value, low_value in zip(high_values, low_values):
                    min_high = min(min_high, high_value)
                    max_low = max(max_low, low_value)
                    if min_high - max_low <= best_score:
                        break
                if min_high - max_low > best_score:
                    best_score = min_high - max_low
                    best_transform = transform
            if best_score > 0:
                transforms.append(best_transform)
        return transforms


def _rotations(step: float) -> Iterator[gemmi.Mat33]:
    for beta in arange(step / 2, 180, step):
        spl = 360 / floor(cos(0.5 * radians(beta)) * 360 / step + 1)
        smi = 360 / floor(sin(0.5 * radians(beta)) * 360 / step + 1)
        for thpl in arange(spl / 2, 720, spl):
            for thmi in arange(smi / 2, 360, smi):
                alpha = 0.5 * (thpl + thmi) % 360
                gamma = 0.5 * (thpl - thmi) % 360
                if alpha <= 360 and beta <= 180 and gamma <= 360:
                    rotation = Rotation.from_euler("ZYZ", [alpha, beta, gamma], True)
                    yield gemmi.Mat33(rotation.as_matrix())
