import gemmi
import numpy as np

from modelcraft.build.sequence import assign, assign_probabilities, residue_numbers
from modelcraft.build.side import LABELS
from modelcraft.tests.ccp4 import insulin_contents, insulin_refmac

CODES = [gemmi.find_tabulated_residue(name).one_letter_code.upper() for name in LABELS]
SEQUENCE = "".join(np.random.default_rng(42).choice(CODES, 300))


def _probabilities(codes, rng, right=0.6):
    "Noisy probabilities with the given probability on the right type"
    noise = rng.dirichlet(np.ones(20), size=len(codes)) * (1 - right)
    for row, code in zip(noise, codes):
        row[CODES.index(code)] += right
    return noise


def test_chain_with_an_insertion_and_a_deletion():
    rng = np.random.default_rng(0)
    # Residues 100-129 with 115 missing and an extra residue after 120
    positions = list(range(100, 115)) + list(range(116, 121)) + [-1]
    positions += list(range(121, 130))
    codes = [SEQUENCE[p] if p >= 0 else "G" for p in positions]
    probabilities = _probabilities(codes, rng)
    probabilities[positions.index(-1)] = rng.dirichlet(np.ones(20))
    (assignment,) = assign_probabilities([probabilities], [SEQUENCE])
    assert np.all(assignment[:, 0][assignment[:, 1] >= 0] == 0)
    expected = np.array(positions)
    assigned = assignment[:, 1] >= 0
    assert assigned.sum() >= 25
    assert np.all(assignment[assigned, 1] == expected[assigned])


def test_chain_without_a_match_is_not_assigned():
    rng = np.random.default_rng(1)
    probabilities = rng.dirichlet(np.ones(20), size=30)
    (assignment,) = assign_probabilities([probabilities], [SEQUENCE])
    assert np.all(assignment == -1)


def test_alignments_do_not_cross_between_sequences():
    rng = np.random.default_rng(2)
    first, second = SEQUENCE[:60], SEQUENCE[60:120]
    codes = first[-10:] + second[:10]
    (assignment,) = assign_probabilities(
        [_probabilities(codes, rng, right=0.9)], [first, second]
    )
    assigned = assignment[assignment[:, 1] >= 0]
    assert len(set(assigned[:, 0])) <= 1


def test_residue_numbers():
    assert residue_numbers([-1, -1]) == [(1, " "), (2, " ")]
    numbers = residue_numbers([-1, 4, -1, 6, -1, -1, -1, 8, -1])
    assert numbers == [
        (4, " "),
        (5, " "),
        (6, " "),
        (7, " "),
        (8, " "),
        (8, "A"),
        (8, "B"),
        (9, " "),
        (10, " "),
    ]


def test_insulin_is_sequenced():
    refmac = insulin_refmac()
    grid = refmac.fphi_best.map(spacing=0.5)
    sequences = [protein.sequence for protein in insulin_contents().proteins]
    chains, numbers = [], []
    for chain in refmac.structure[0]:
        residues = [
            r for r in chain if all(r.find_atom(a, "*") for a in ("N", "CA", "C"))
        ]
        chains.append(
            np.array(
                [[r.find_atom(a, "*").pos.tolist() for a in ("N", "CA", "C")]
                 for r in residues]
            )
        )  # fmt: skip
        numbers.append(np.array([r.seqid.num for r in residues]))
    right = total = 0
    for c, (assignment, nums) in enumerate(
        zip(assign(grid, chains, sequences), numbers)
    ):
        assert np.all(assignment[assignment[:, 1] >= 0, 0] == c)
        right += np.sum(assignment[:, 1] + 1 == nums)
        total += len(nums)
    assert right / total > 0.8
