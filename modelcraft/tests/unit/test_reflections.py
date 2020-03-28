from modelcraft.tests import data_path
from modelcraft.reflections import DataFile


def test_1kv9():
    hklin = DataFile(data_path("1kv9_data.mtz"))
    assert round(hklin.cell.a, 3) == 54.807
    assert round(hklin.cell.b, 3) == 57.437
    assert round(hklin.cell.c, 3) == 67.523
    assert round(hklin.cell.alpha, 2) == 89.65
    assert round(hklin.cell.beta, 2) == 69.34
    assert round(hklin.cell.gamma, 2) == 68.39
    assert hklin.spacegroup.hm == "P 1"
    assert hklin.num_reflections == 51422
    assert round(hklin.resolution_high, 2) == 1.80
    assert round(hklin.resolution_low, 2) == 29.87
    columns = ["H", "K", "L", "FREE", "FP", "SIGFP", "HLA", "HLB", "HLC", "HLD"]
    assert len(hklin.columns) == len(columns)
    assert all(col in hklin for col in columns)
    assert ",".join(columns) in hklin
    assert "LABEL_NOT_IN_FILE" not in hklin
    assert hklin.frees == {"FREE"}
    assert hklin.fsigfs == {"FP,SIGFP"}
    assert hklin.abcds == {"HLA,HLB,HLC,HLD"}
    assert hklin.phifoms == set()
    assert hklin.fphis == set()


def test_hewl():
    hklin = DataFile(data_path("hewl_data.mtz"))
    assert hklin.spacegroup.hm == "P 43 21 2"
    assert hklin.num_reflections == 10055
    assert round(hklin.resolution_high, 2) == 1.86
    assert round(hklin.resolution_low, 2) == 55.23
    columns = [
        "H",
        "K",
        "L",
        "F_New",
        "SIGF_New",
        "DANO_New",
        "SIGDANO_New",
        "F_New(+)",
        "SIGF_New(+)",
        "F_New(-)",
        "SIGF_New(-)",
        "IMEAN_New",
        "SIGIMEAN_New",
        "I_New(+)",
        "SIGI_New(+)",
        "I_New(-)",
        "SIGI_New(-)",
        "ISYM_New",
        "FreeR_flag",
        "FWT",
        "PHWT",
        "PHIB",
        "FOM",
        "FPFOM",
        "HLA",
        "HLB",
        "HLC",
        "HLD",
        "HLanomA",
        "HLanomB",
        "HLanomC",
        "HLanomD",
        "parrot.ABCD.A",
        "parrot.ABCD.B",
        "parrot.ABCD.C",
        "parrot.ABCD.D",
        "parrot.F_phi.F",
        "parrot.F_phi.phi",
    ]
    assert len(hklin.columns) == len(columns)
    assert all(col in hklin for col in columns)
    assert hklin.frees == {"FreeR_flag"}
    assert hklin.fsigfs == {"F_New,SIGF_New"}
    assert hklin.abcds == {
        "HLA,HLB,HLC,HLD",
        "HLanomA,HLanomB,HLanomC,HLanomD",
        "parrot.ABCD.A,parrot.ABCD.B,parrot.ABCD.C,parrot.ABCD.D",
    }
    assert hklin.phifoms == {
        "PHWT,FPFOM",
        "parrot.F_phi.phi,FPFOM",
        "PHIB,FPFOM",
        "PHIB,FOM",
        "parrot.F_phi.phi,FOM",
        "PHWT,FOM",
    }
    assert hklin.fphis == {
        "FWT,PHWT",
        "parrot.F_phi.F,parrot.F_phi.phi",
    }
