import pytest
from modelcraft.mcmodel import _slice_min


@pytest.mark.parametrize(
    "chain_rsccs, fragment_rsccs, start_index, expected",
    [
        ([0.11, 0.11, 0.11], [0.99, 0.99], 1, 0),
        ([0.11, 0.11, None], [0.99, 0.99], 1, 0),
        ([0.11, None, None], [0.99, 0.99], 1, 0),
        ([None, None, None], [0.99, 0.99], 1, 0),
        ([None, None, 0.11], [0.99, 0.99], 1, 0),
        ([None, 0.11, 0.11], [0.99, 0.99], 1, 0),
        ([0.99, 0.11, 0.11], [0.99, 0.99], 1, 0),
        ([0.99, 0.99, 0.11], [0.11, 0.99], 1, 1),
        ([0.99, 0.99, 0.99], [0.11, 0.11], 1, 2),
        ([0.11, 0.11, 0.11], [0.99, 0.99, 0.99], 0, 0),
        ([0.11, 0.99, 0.99], [0.99, 0.11, 0.11], 0, 3),
    ],
)
def test_slice_min(chain_rsccs, fragment_rsccs, start_index, expected):
    assert _slice_min(chain_rsccs, fragment_rsccs, start_index) == expected
