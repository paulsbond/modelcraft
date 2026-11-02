import numpy as np
import pytest

from modelcraft.build.grow import next_residues, prev_residues
from modelcraft.tests.ccp4 import insulin_refmac

NAMES = ("N", "CA", "C")


def _chains(reverse=False):
    "Insulin's chains as (n, 3, 3) arrays of N, CA and C, and the map"
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    chains = [
        np.array(
            [
                [r.find_atom(a, "*").pos.tolist() for a in NAMES]
                for r in chain
                if all(r.find_atom(a, "*") for a in NAMES)
            ]
        )
        for chain in refmac.structure[0]
    ]
    return [c[::-1] for c in chains] if reverse else chains, grid


def _max_errors(placed, target) -> np.ndarray:
    return np.linalg.norm(placed - target, axis=-1).max(axis=-1)


@pytest.mark.parametrize("grow", [next_residues, prev_residues])
def test_grown_residues_land_on_the_model(grow):
    chains, grid = _chains(reverse=grow is prev_residues)
    right = total = 0
    for chain in chains:
        placed, errors = grow(grid, chain[:-1])
        assert np.all((errors > 0) & (errors < 2))
        right += np.sum(_max_errors(placed, chain[1:]) < 1)
        total += len(placed)
    assert total > 40
    assert right / total > 0.75


@pytest.mark.parametrize("grow", [next_residues, prev_residues])
def test_grown_chains_follow_the_model(grow):
    chains, grid = _chains(reverse=grow is prev_residues)
    fractions = []
    for chain in chains:
        for i in range(len(chain) - 1):
            current, steps = chain[i : i + 1], 0
            for target in chain[i + 1 :]:
                current, _ = grow(grid, current)
                if _max_errors(current[0], target) > 1.5:
                    break
                steps += 1
            fractions.append(steps / (len(chain) - 1 - i))
    assert np.mean(fractions) > 0.35


def test_no_residues():
    _, grid = _chains()
    placed, errors = next_residues(grid, np.zeros((0, 3, 3)))
    assert placed.shape == (0, 3, 3) and errors.shape == (0,)
