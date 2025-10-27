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

## Graphical interfaces

For X-ray crystallography, ModelCraft can be used through
[CCP4i2](https://ccp4i2.gitlab.io/rstdocs/tasks/modelcraft/index.html) or
[CCP4 Cloud](https://cloud.ccp4.ac.uk/manuals/html-taskref/doc.task.ModelCraft.html).
For cryo-EM, it is available in
[Doppio](https://www.ccpem.ac.uk/docs/doppio/user_guide.html).

## Command-line interface

The first argument must be either
`xray` for X-ray crystallography or `em` for cryo-EM SPA.
The command line documentation has detailed information on individual arguments.

```bash
modelcraft xray --help
modelcraft em --help
```

## Links

- [Source (latest)](https://github.com/paulsbond/modelcraft)
- [Source (releases)](https://github.com/paulsbond/modelcraft/releases)
- [Issues (e.g. bugs and feature requests)](https://github.com/paulsbond/modelcraft/issues)
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
