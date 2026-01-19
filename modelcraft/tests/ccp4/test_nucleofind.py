import gemmi

from modelcraft.jobs.freerflag import FreeRFlag
from modelcraft.jobs.nucleofind import NucleoFind
from modelcraft.jobs.refmac import Refmac
from modelcraft.pipeline import Pipeline
from modelcraft.reflections import DataItem
from modelcraft.scripts.contents import _entry_contents
from modelcraft.structure import (
    ModelStats,
    contains_residue,
    read_structure,
    remove_residues,
)

from ...jobs.nucleofind_build import NucleoFindBuild
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
    refmac = Refmac(structure=structure, fsigf=fsigf, freer=freer, cycles=0).run()
    contents = _entry_contents("102d")
    pipeline = Pipeline(keep_jobs=True)
    # Test without an input structure
    nucleofind_result = NucleoFind(
        fphi=refmac.fphi_best,
    ).run(pipeline)

    build_result = NucleoFindBuild(
        contents=contents,
        fsigf=fsigf,
        phases=refmac.abcd,
        fphi=refmac.fphi_best,
        freer=freer,
        nucleofind_result=nucleofind_result,
    ).run(pipeline)

    stats = ModelStats(build_result.structure)
    assert stats.residues > 12

    # Test with an input structure
    remove_residues(structure, ["HOH"])
    nucleofind_result = NucleoFind(
        fphi=refmac.fphi_best,
    ).run(pipeline)

    build_result = NucleoFindBuild(
        contents=contents,
        fsigf=fsigf,
        phases=refmac.abcd,
        fphi=refmac.fphi_best,
        freer=freer,
        structure=structure,
        nucleofind_result=nucleofind_result,
    ).run(pipeline)
    stats = ModelStats(build_result.structure)
    assert contains_residue(build_result.structure, "DT")
    assert stats.residues > 22
    assert build_result.fragments_built > 2
    assert build_result.residues_built > 22
    assert build_result.residues_sequenced > 20
    assert build_result.longest_fragment > 11
    # TODO: assert contains_residue(nautilus.structure, "TNT")
