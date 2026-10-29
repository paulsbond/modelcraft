"""Degrade reference data to resemble work data (Buccaneer's MapSimulate).

The reference amplitudes are rescaled to the work data's fall-off with resolution.
Within resolution bins the work FOMs are transferred by rank of |E| onto the
reference reflections, and phase errors are drawn to match those FOMs.
"""

import dataclasses

import gemmi
import numpy as np
from scipy.special import i0e, i1e

from .reflections import DataItem

_KAPPA = np.concatenate([[0.0], np.logspace(-4, 4, 2000)])
_SIM = i1e(_KAPPA) / i0e(_KAPPA)


def _invsim(fom: np.ndarray) -> np.ndarray:
    "Von Mises concentration with I1/I0 equal to each FOM"
    return np.interp(fom, _SIM, _KAPPA)


def phases_from_hl(hl: np.ndarray, centric: np.ndarray):
    "Centroid phases (radians) and FOMs from Hendrickson-Lattman coefficients"
    phi = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    cos1, sin1, cos2, sin2 = np.cos(phi), np.sin(phi), np.cos(2 * phi), np.sin(2 * phi)
    exponent = np.outer(hl[:, 0], cos1) + np.outer(hl[:, 1], sin1)
    exponent += np.outer(hl[:, 2], cos2) + np.outer(hl[:, 3], sin2)
    exponent -= exponent.max(axis=1, keepdims=True)
    probability = np.exp(exponent)
    centroid = probability @ np.exp(1j * phi) / probability.sum(axis=1)
    fom = np.abs(centroid)
    fom[centric] = np.tanh(np.hypot(hl[centric, 0], hl[centric, 1]))
    return np.angle(centroid), np.nan_to_num(fom)


def _arrays(item: DataItem):
    data = np.array(item, copy=True)
    return data[:, :3].astype(int), data[:, 3:]


def _phases_and_foms(phases: DataItem):
    "Phases (radians) and FOMs from HL or PHI,FOM columns"
    hkl, values = _arrays(phases)
    if phases.types == "AAAA":
        centric = phases.spacegroup.operations().centric_flag_array(hkl)
        return hkl, *phases_from_hl(values, centric.astype(bool))
    return hkl, np.radians(values[:, 0]), np.nan_to_num(values[:, 1])


@dataclasses.dataclass
class _Data:
    "Observed reflections to a resolution with their phases and FOMs"

    hkl: np.ndarray
    f: np.ndarray
    s: np.ndarray
    centric: np.ndarray
    esq: np.ndarray
    phi: np.ndarray
    fom: np.ndarray

    @classmethod
    def from_items(cls, fsigf: DataItem, phases: DataItem, resolution: float):
        hkl, values = _arrays(fsigf)
        f = values[:, 0]
        s = fsigf.cell.calculate_1_d2_array(hkl)
        keep = (s <= 1 / resolution**2) & ~np.isnan(f)
        hkl, f, s = hkl[keep], f[keep], s[keep]
        ops = fsigf.spacegroup.operations()
        phase_hkl, phi, fom = _phases_and_foms(phases)
        lookup = {tuple(h): i for i, h in enumerate(phase_hkl)}
        index = np.array([lookup.get(tuple(h), -1) for h in hkl])
        return cls(
            hkl=hkl,
            f=f,
            s=s,
            centric=ops.centric_flag_array(hkl).astype(bool),
            esq=f**2 / ops.epsilon_factor_array(hkl),
            phi=np.where(index >= 0, phi[index], np.nan),
            fom=np.where(index >= 0, np.nan_to_num(fom[index]), 0.0),
        )

    def mean_esq(self, nbins: int = 20):
        "Mean |E|^2 as a function of 1/d^2 (linear between equal-count bins)"
        bins = np.array_split(np.argsort(self.s), nbins)
        centres = np.array([self.s[b].mean() for b in bins])
        means = np.array([self.esq[b].mean() for b in bins])
        return lambda x: np.interp(x, centres, means)

    def bins(self, nbins: int):
        "Acentric and centric indices in equal-count resolution bins"
        ranks = np.empty(len(self.s))
        ranks[np.argsort(self.s, kind="stable")] = np.arange(len(self.s)) / len(self.s)
        index = np.minimum((ranks * nbins).astype(int), nbins - 1)
        acen = [np.flatnonzero((index == i) & ~self.centric) for i in range(nbins)]
        ccen = [np.flatnonzero((index == i) & self.centric) for i in range(nbins)]
        return acen, ccen


def _pad(bins, bin_min: int):
    "Bins extended with alternately higher and lower neighbours to a minimum size"
    padded = []
    for i in range(len(bins)):
        indices, j = list(bins[i]), 1
        while len(indices) < bin_min and j < len(bins):
            if 0 <= i + j < len(bins):
                indices.extend(bins[i + j])
            j = -j if j > 0 else -j + 1
        padded.append(np.array(indices, dtype=int))
    return padded


def _transfer_foms(ref: _Data, work: _Data, nbins: int, bin_min: int):
    "Work FOMs given to the reference reflections by rank of |E| in each bin"
    r_acen, r_ccen = ref.bins(nbins)
    w_acen, w_ccen = work.bins(nbins)
    w_ccen = [c if len(c) > 0 else a for c, a in zip(_pad(w_ccen, bin_min), w_acen)]
    r_ccen = [c if len(c) > 0 else a for c, a in zip(r_ccen, r_acen)]
    fom = ref.fom.copy()
    valid = np.ones(len(ref.f), dtype=bool)
    for r_bins, w_bins in ((r_acen, w_acen), (r_ccen, w_ccen)):
        for r_bin, w_bin in zip(r_bins, w_bins):
            if len(w_bin) == 0:
                valid[r_bin] = False
            elif len(r_bin) > 0:
                r_sorted = r_bin[np.argsort(ref.esq[r_bin], kind="stable")]
                w_sorted = w_bin[np.argsort(work.esq[w_bin], kind="stable")]
                k = np.arange(len(r_sorted)) * len(w_sorted) // len(r_sorted)
                fom[r_sorted] = work.fom[w_sorted[k]]
    return fom, valid


def simulate(
    ref_fsigf: DataItem,
    ref_phases: DataItem,
    work_fsigf: DataItem,
    work_phases: DataItem,
    resolution: float,
    *,
    nbins: int = 100,
    bin_min: int = 20,
    seed: int = 0,
) -> gemmi.FloatGrid:
    "Reference map with the resolution, fall-off and phase errors of the work data"
    ref = _Data.from_items(ref_fsigf, ref_phases, resolution)
    work = _Data.from_items(work_fsigf, work_phases, resolution)

    # Scale the reference amplitudes to the fall-off of the work data
    cellscale = ref_fsigf.cell.volume / work_fsigf.cell.volume
    f = ref.f * np.sqrt(cellscale * work.mean_esq()(ref.s) / ref.mean_esq()(ref.s))
    fom, valid = _transfer_foms(ref, work, nbins, bin_min)

    # Phase errors to go with the FOMs
    rng = np.random.default_rng(seed)
    dphi = rng.vonmises(0.0, np.maximum(_invsim(fom), 1e-8))
    flip = rng.random(len(f)) >= 0.5 * (1 + fom)
    dphi[ref.centric] = np.where(flip[ref.centric], np.pi, 0.0)

    keep = valid & ~np.isnan(ref.phi)
    mtz = gemmi.Mtz(with_base=True)
    mtz.cell = ref_fsigf.cell
    mtz.spacegroup = ref_fsigf.spacegroup
    mtz.add_column("F", "F")
    mtz.add_column("PHI", "P")
    phi = np.degrees(ref.phi + dphi)
    mtz.set_data(np.column_stack([ref.hkl[keep], (f * fom)[keep], phi[keep]]))
    return mtz.transform_f_phi_to_map("F", "PHI", sample_rate=3.0)
