import gemmi
import numpy as np
import scipy.linalg


def remove_scale(structure: gemmi.Structure) -> None:
    "Convert the structure to use standard conventions for axes"
    new_cell = gemmi.UnitCell(*structure.cell.parameters)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    fractional = structure.cell.fractionalize(atom.pos)
                    atom.pos = new_cell.orthogonalize(fractional)
    structure.cell = new_cell


def max_distortion(old_cell: gemmi.UnitCell, new_cell: gemmi.UnitCell) -> float:
    "Return the maximum distortion between two cells as a fraction"
    frac1 = old_cell.frac.mat
    orth2 = new_cell.orth.mat
    identity = gemmi.Mat33()
    matrix = np.array(identity) - np.array(orth2.multiply(frac1))
    eigenvalues, _ = scipy.linalg.eig(matrix)
    return max(abs(value) for value in eigenvalues.real)


def update_cell(structure: gemmi.Structure, new_cell: gemmi.UnitCell) -> None:
    "Update the structure cell without distorting the model structure"
    for model in structure:
        atoms = [atom for chain in model for residue in chain for atom in residue]
        if len(atoms) > 2:
            old_positions = []
            aim_positions = []
            for atom in atoms:
                old_positions.append(atom.pos)
                fractional = structure.cell.fractionalize(atom.pos)
                orthogonal = new_cell.orthogonalize(fractional)
                aim_positions.append(orthogonal)
            result = gemmi.superpose_positions(aim_positions, old_positions)
            if not np.isnan(result.rmsd):
                for atom in atoms:
                    atom.pos = gemmi.Position(result.transform.apply(atom.pos))
        else:
            for atom in atoms:
                fractional = structure.cell.fractionalize(atom.pos)
                orthogonal = new_cell.orthogonalize(fractional)
                atom.pos = orthogonal
    structure.cell = new_cell
