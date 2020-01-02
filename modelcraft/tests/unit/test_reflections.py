from modelcraft.tests import data_path
from modelcraft.reflections import ReflectionFile


def test_resolution():
    hklin = ReflectionFile(data_path("1kv9_data.mtz"))
    assert hklin.resolution() == 1.80
