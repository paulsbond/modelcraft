"""Assign sequence to built chains from the side network's residue type probabilities.

Each chain is aligned to the sequences with a local profile HMM that allows for
insertions and deletions in the model.
"""

import gemmi
import numpy as np
from scipy.special import logsumexp

from . import side
from .side import LABELS

_ONE_LETTER = [
    gemmi.find_tabulated_residue(name).one_letter_code.upper() for name in LABELS
]


def _sequence_arrays(sequences):
    "Type index in LABELS (-1 if unknown), sequence index and position of each residue"
    types, which, position = [], [], []
    for s, sequence in enumerate(sequences):
        for j, code in enumerate(sequence.upper()):
            types.append(_ONE_LETTER.index(code) if code in _ONE_LETTER else -1)
            which.append(s)
            position.append(j)
    return np.array(types, dtype=int), np.array(which), np.array(position)


def _log_odds(probabilities, types, composition, mix):
    """Log-odds (n, m) of built residues being each sequence residue.

    The probabilities are mixed with the sequence composition so that a confidently
    wrong residue costs less than an insertion and deletion around it.
    """
    mixed = (1 - mix) * probabilities + mix * composition
    scores = np.log(mixed) - np.log(composition)
    result = scores[:, np.clip(types, 0, None)]
    result[:, types < 0] = 0.0
    return result


def _align(
    scores,
    which,
    *,
    start=0.01,
    end=0.01,
    insertion=0.05,
    extend=0.5,
    deletion=0.02,
    max_skip=3,
):
    """Posterior probabilities (n, m) of each built residue being at each position.

    which (m,) gives the sequence of each position, so that alignments stay within
    one sequence. Deletions skip up to max_skip positions.
    """
    n, m = scores.shape
    ninf = -np.inf
    log_match = np.log(1 - insertion - max_skip * deletion - end)
    log_start, log_stay = np.log(start / m), np.log(1 - start)
    log_ins, log_ext, log_resume = np.log(insertion), np.log(extend), np.log(1 - extend)
    log_del, log_end = np.log(deletion), np.log(end)
    # Moves from position j to j + 1 + k within one sequence, as (m,) masks on j
    same = []
    for k in range(max_skip + 1):
        mask = np.zeros(m, bool)
        if 1 + k < m:
            mask[: m - 1 - k] = which[1 + k :] == which[: m - 1 - k]
        same.append(mask)
    moves = [
        (k, same[k], log_match if k == 0 else log_del) for k in range(max_skip + 1)
    ]

    def forward_shift(values, k):
        "values[j - k] at j"
        out = np.full(m, ninf)
        if k < m:
            out[k:] = values[: m - k]
        return out

    def backward_shift(values, k):
        "values[j + k] at j"
        out = np.full(m, ninf)
        if k < m:
            out[: m - k] = values[k:]
        return out

    # Forward: match (fm), insertion (fi), before (fb) and after (fa) the alignment
    fm, fi = np.full((n, m), ninf), np.full((n, m), ninf)
    fb, fa = np.full(n, ninf), np.full(n, ninf)
    fb[0] = log_stay
    fm[0] = log_start + scores[0]
    for r in range(1, n):
        fb[r] = fb[r - 1] + log_stay
        terms = [np.full(m, fb[r - 1] + log_start)]
        for k, mask, log_move in moves:
            terms.append(
                forward_shift(np.where(mask, fm[r - 1], ninf), 1 + k) + log_move
            )
        terms.append(forward_shift(np.where(same[0], fi[r - 1], ninf), 1) + log_resume)
        fm[r] = logsumexp(np.stack(terms), axis=0) + scores[r]
        fi[r] = np.logaddexp(fm[r - 1] + log_ins, fi[r - 1] + log_ext)
        fa[r] = np.logaddexp(fa[r - 1], logsumexp(fm[r - 1]) + log_end)
    total = logsumexp([fb[-1], fa[-1], logsumexp(fm[-1]), logsumexp(fi[-1])])

    # Backward (after the alignment has a backward value of 0)
    bm, bi = np.zeros((n, m)), np.zeros((n, m))
    for r in range(n - 2, -1, -1):
        following = bm[r + 1] + scores[r + 1]
        terms = [
            np.where(mask, backward_shift(following, 1 + k), ninf) + log_move
            for k, mask, log_move in moves
        ]
        terms.append(bi[r + 1] + log_ins)
        terms.append(np.full(m, log_end))
        bm[r] = logsumexp(np.stack(terms), axis=0)
        resume = np.where(same[0], backward_shift(following, 1), ninf) + log_resume
        bi[r] = np.logaddexp(bi[r + 1] + log_ext, resume)
    return np.exp(fm + bm - total)


def assign(grid, chains, sequences, **options):
    """Sequence (index, position) arrays (n, 2) for each chain, -1 where unassigned.

    The options are as for assign_probabilities.
    """
    probabilities = [
        side.type_probabilities(side.probabilities(grid, chain)) for chain in chains
    ]
    return assign_probabilities(probabilities, sequences, **options)


def assign_probabilities(
    probabilities, sequences, mix=0.3, min_probability=0.5, **options
):
    """As assign, from type probabilities (n, 20) for each chain.

    A residue is assigned its most probable position if that probability is above
    min_probability and the position follows those of the residues before it.
    """
    types, which, position = _sequence_arrays(sequences)
    composition = np.bincount(types[types >= 0], minlength=len(LABELS)) + 1.0
    composition /= composition.sum()
    result = []
    for chain in probabilities:
        assignment = np.full((len(chain), 2), -1)
        if len(chain) and len(types):
            scores = _log_odds(chain, types, composition, mix)
            posterior = _align(scores, which, **options)
            best = posterior.argmax(axis=1)
            last = -1
            for r, j in enumerate(best):
                if posterior[r, j] > min_probability and position[j] > last:
                    assignment[r] = which[j], position[j]
                    last = position[j]
        result.append(assignment)
    return result


def residue_numbers(positions):
    """Residue numbers and insertion codes for a chain from its sequence positions
    (-1 where unassigned), with insertion codes where there aren't enough numbers."""
    n = len(positions)
    assigned = [r for r in range(n) if positions[r] >= 0]
    if not assigned:
        return [(r + 1, " ") for r in range(n)]
    numbers = [None] * n
    for r in assigned:
        numbers[r] = (int(positions[r]) + 1, " ")
    first, last = assigned[0], assigned[-1]
    for r in range(first):
        numbers[r] = (numbers[first][0] - (first - r), " ")
    for r in range(last + 1, n):
        numbers[r] = (numbers[last][0] + (r - last), " ")
    for a, b in zip(assigned, assigned[1:]):
        low, high = numbers[a][0], numbers[b][0]
        codes = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        for r in range(a + 1, b):
            if low + (r - a) < high:
                numbers[r] = (low + (r - a), " ")
            else:
                numbers[r] = (high - 1, next(codes, " "))
    return numbers
