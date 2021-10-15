import gemmi
from modelcraft.jobs.freerflag import FreeRFlag
from modelcraft.jobs.nautilus import Nautilus
from modelcraft.jobs.refmac import RefmacXray
from modelcraft.pipeline import Pipeline
from modelcraft.reflections import DataItem
from modelcraft.scripts.contents import _entry_contents
from modelcraft.structure import (
    contains_residue,
    ModelStats,
    read_structure,
    remove_residues,
)
from . import in_temp_directory, pdbe_download


@in_temp_directory
def test_102d():
    # Prepare input data
    pdbe_download("102d.cif")
    structure = read_structure("102d.cif")
    pdbe_download("r102dsf.ent")
    doc = gemmi.cif.read("r102dsf.ent")
    rblocks = gemmi.as_refln_blocks(doc)
    cif2mtz = gemmi.CifToMtz()
    mtz = cif2mtz.convert_block_to_mtz(rblocks[0])
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = FreeRFlag(fsigf).run().freer
    refmac = RefmacXray(structure=structure, fsigf=fsigf, freer=freer, cycles=0).run()
    contents = _entry_contents("102d")
    pipeline = Pipeline(keep_jobs=True)
    # Test without an input structure
    nautilus = Nautilus(
        contents=contents,
        fsigf=fsigf,
        phases=refmac.abcd,
        fphi=refmac.fphi_best,
        freer=freer,
    ).run(pipeline)
    stats = ModelStats(nautilus.structure)
    assert stats.residues > 12
    # Test with an input structure
    remove_residues(structure, ["HOH"])
    nautilus = Nautilus(
        contents=contents,
        fsigf=fsigf,
        phases=refmac.abcd,
        fphi=refmac.fphi_best,
        freer=freer,
        structure=structure,
    ).run(pipeline)
    stats = ModelStats(nautilus.structure)
    assert contains_residue(nautilus.structure, " DT")
    assert stats.residues > 22
    # TODO: assert contains_residue(nautilus.structure, "TNT")
