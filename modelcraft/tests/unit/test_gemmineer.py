from modelcraft.tests import data_path
import modelcraft.coordinates as coordinates


def test_number_of_known():
    assert len(coordinates._known_protein_residues) == 22


def test_unk_is_known():
    assert "UNK" in coordinates._known_protein_residues


def test_1kv9_model_stats():
    xyzin = data_path("1kv9_model.pdb")
    coordinates.CoordinateFile(xyzin)
