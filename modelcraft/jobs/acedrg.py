from typing import Optional
import gemmi
from .job import Job


class Acedrg(Job):
    def __init__(
        self,
        code: str,
        cif: Optional[gemmi.cif.Document] = None,
        smiles: Optional[str] = None,
    ):
        super().__init__("acedrg")

        code = code.upper()
        args = ["-r", code]

        if cif is not None:
            cif.write_file(self.path("input.cif"))
            args += ["-c", self.path("input.cif")]

        if smiles is not None:
            args += ["-i", smiles]

        args += ["-o", self.path("acedrg")]
        self.run("acedrg", args)

        acedrg_cif = self.path("acedrg.cif")
        block = gemmi.cif.read(acedrg_cif)[-1]
        self.chemcomp = gemmi.make_chemcomp_from_block(block)

        self.finish()
