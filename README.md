![](https://github.com/paulsbond/modelcraft/blob/691b8a0d72dbf48b8d8a431bdb193298a397e1bb/modelcraft.png)

# ModelCraft

ModelCraft is an automated model-building pipeline
that builds proteins and nucleic acids
into X-ray crystallography and cryo-EM maps.
It combines the previous Buccaneer and Nautilus pipelines and adds new steps
for density modification, refinement, validation and pruning.

## Installation

ModelCraft is distributed with CCP4.
Please obtain the latest version from the
[CCP4 download page](https://www.ccp4.ac.uk/download).
To build into cryo-EM maps, you must also
[install CCP-EM](https://www.ccpem.ac.uk/download.php)
and source both the CCP4 and CCP-EM environments.

## Graphical interfaces

For X-ray crystallography, ModelCraft can be used through
[CCP4i2](https://ccp4i2.gitlab.io/rstdocs/tasks/modelcraft/index.html) or
[CCP4 Cloud](https://cloud.ccp4.ac.uk/manuals/html-taskref/doc.task.ModelCraft.html).
For cryo-EM, it is available in
[Doppio](https://www.ccpem.ac.uk/docs/doppio/user_guide.html).

## Command-line interface

The first argument must be either
`xray` for X-ray crystallography or `em` for cryo-EM.

The simplest execution for X-ray crystallography requires only a
[description of the asymmetric unit contents](#asu-contents-description)
and a reflection data file in MTZ format
with observations, a free-R flag and starting phases,
e.g. from experimental phasing.

```bash
modelcraft xray --contents contents.json --data data.mtz
```

More commonly, you will have an partial atomic model that you want to complete,
e.g. a placed molecular replacement model.
A starting model can be provided in either PDB, mmCIF or mmJSON format.
It will be refined and used as a starting point
instead of starting from phases in the data file.

```bash
modelcraft xray --contents contents.json --data data.mtz --model model.cif
```

For cryo-EM, either two halfmaps or a single map must be provided,
along with a resolution.
A starting model can be provided in the same way as with X-ray data.

```bash
modelcraft em --contents contents.json --map half1.mrc half2.mrc --resolution 2.5
modelcraft em --contents contents.json --map map.mrc --resolution 2.5
```

The command line documentation
has more detailed information on individual arguments.

```bash
modelcraft xray --help
modelcraft em --help
```

### ASU contents description

A description of the expected contents of the asymmetric unit
must be provided as either a FASTA sequence file or a JSON file
using the `--contents` argument.
A sequence file is simpler,
but the JSON format has the following advantages:

- Number of copies and stoichiometry can be specified
  for a more accurate calculation of the solvent fraction.
- Carbohydrates, ligands and buffers may be specified
  in addition to protein, RNA and DNA.
- Molecule types do not need to be guessed from
  sequences (which may be ambiguous).

In order to create a JSON file it may be helpful
to start from the contents for an existing PDB entry.
The `modelcraft-contents` script
creates a contents JSON file for a released PDB entry.

An example JSON file is shown below:

```json
{
"copies": 2,
"proteins": [
    {
        "sequence": "LPGECSVNVIPKMNLDKAKFFSGTWYETHYLDMDPQATEKFCFSFAPRESGGTVMEALYHFNVDSKV",
        "stoichiometry": 1,
        "modifications": ["M->MSE"]
    },
    {
        "sequence": "GGG"
    }
],
"rnas": [
    {
        "sequence": "GGUAACUGUUACAGUUACC",
        "stoichiometry": 2,
        "modifications": ["1->GTP", "19->CCC"]
    }
],
"dnas": [],
"carbs": [
    { "codes": { "NAG": 2 }, "stoichiometry": 1 },
    { "codes": { "MAN": 1, "NAG": 2 }, "stoichiometry": 1 }
],
"ligands": [
    { "code": "HEM", "stoichiometry": 1 }
],
"buffers": ["GOL", "NA", "CL"]
}
```

The file has a list of `proteins`, `rnas`, `dnas`, `carbs`, `ligands`,
and `buffers` that are in the crystal.
The only mandatory items are that each
protein, RNA or DNA chain must have a `sequence`,
each carbohydrate must have a dictionary of `codes`
to specify the number of each sugar,
and each ligand must have a single `code`.

Each component (other than buffers) has a
`stoichiometry` parameter to specify the stoichiometry.
In the example above there are 2 RNA chains for each protein chain.
If the stoichiometry is not specified it is assumed to be 1.
There is also a `copies` parameter for the whole file
to specify how many copies of the contents are in the asymmetric unit.
If this value is not known the most likely number will be estimated.
The `modelcraft-copies` script can be used to
view the solvent fraction and probability for each number of copies
given a contents file and an MTZ file.
It is assumed that the number of ordered buffer molecules is unknown
so they are not included in the solvent calculation.

Finally, protein, RNA and DNA chains may have
a list of `modifications`,
e.g. `M->MSE` to specify that
all methionine residues are actually selenomethionine
or `1->GTP` to specify that
the residue 1 is guanosine triphosphate.

**Note:**
ModelCraft does not yet build carbohydrates, ligands,
or modified residues (other than selenomethionine derivatives).
However, this is planned for the future
and inclusion of these components in the contents
allows for more accurate calculation of the solvent fraction
during Parrot density modification in the X-ray pipeline.

## Links

- [Source (latest)](https://github.com/paulsbond/modelcraft)
- [Source (releases)](https://github.com/paulsbond/modelcraft/releases)
- [Issues (e.g. bugs and feature requests)](https://github.com/paulsbond/modelcraft/issues">)
- [PyPI (Python Package Index)](https://pypi.org/project/modelcraft)

## Citations

- [ModelCraft](https://doi.org/10.1107/S2059798322007732)
- [Buccaneer](https://doi.org/10.1107/S0907444906022116)
- [EMDA](https://doi.org/10.1016/j.jsb.2021.107826)
- [Nautilus](https://doi.org/10.1107/S2052252514019290)
- [Parrot](https://doi.org/10.1107/S090744490903947X)
- [Refmac](https://doi.org/10.1107/S2059798318000979)
- [Servalcat](https://doi.org/10.1107/S2059798321009475)
- [Sheetbend](https://doi.org/10.1107/S2059798320013170)
