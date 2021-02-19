from typing import List, Optional
import json
import modelcraft.pdbe as pdbe
from .polymer import Polymer, PolymerType
from .carb import Carb
from .ligand import Ligand


class AsuContents:
    def __init__(self, path_or_pdbid: Optional[str] = None):
        self.polymers: List[Polymer] = []
        self.carbs: List[Carb] = []
        self.ligands: List[Ligand] = []
        if path_or_pdbid is not None:
            if len(path_or_pdbid) == 4:
                self.add_from_pdbid(path_or_pdbid)
            else:
                if path_or_pdbid[-5:] == ".json":
                    self.add_from_json_file(path_or_pdbid)
                else:
                    self.add_from_sequence_file(path_or_pdbid)

    def add_from_sequence_file(
        self, path: str, polymer_type: Optional[PolymerType] = None
    ) -> None:
        for polymer in Polymer.from_sequence_file(path, polymer_type):
            self.polymers.append(polymer)

    def add_from_json_file(self, path: str) -> None:
        with open(path) as stream:
            components = json.load(stream)
        for component in components:
            if "sequence" in component:
                polymer = Polymer.from_component_json(component)
                self.polymers.append(polymer)
            elif "codes" in component:
                carb = Carb.from_component_json(component)
                self.carbs.append(carb)
            elif "code" in component:
                ligand = Ligand.from_component_json(component)
                self.ligands.append(ligand)

    def add_from_pdbid(self, pdbid: str) -> None:
        for mol in pdbe.molecules(pdbid):
            if "sequence" in mol:
                polymer = Polymer.from_pdbe_molecule_dict(mol)
                self.polymers.append(polymer)
            if mol["molecule_type"] == "carbohydrate polymer":
                carb = Carb.from_pdbe_molecule_dict(mol)
                self.carbs.append(carb)
            if mol["molecule_type"] == "bound":
                ligand = Ligand.from_pdbe_molecule_dict(mol)
                if ligand.code not in ("UNL", "UNX"):
                    self.ligands.append(ligand)

    def is_selenomet(self) -> bool:
        proteins = [p for p in self.polymers if p.type == PolymerType.PROTEIN]
        return len(proteins) > 0 and all(p.is_selenomet() for p in proteins)

    def components_json(self) -> list:
        polymers = [polymer.to_component_json() for polymer in self.polymers]
        carbs = [carb.to_component_json() for carb in self.carbs]
        ligands = [ligand.to_component_json() for ligand in self.ligands]
        return polymers + carbs + ligands

    def write_sequence_file(
        self,
        path: str,
        polymer_type: Optional[PolymerType] = None,
        line_length: int = 60,
    ) -> None:
        with open(path, "w") as stream:
            for polymer in self.polymers:
                if polymer_type is None or polymer.type == polymer_type:
                    stream.write(f">{polymer.type.value}\n")
                    for i in range(0, len(polymer.sequence), line_length):
                        stream.write(polymer.sequence[i : i + line_length] + "\n")

    def write_json_file(self, path: str) -> None:
        with open(path, "w") as stream:
            json.dump(self.components_json(), stream, indent=2)
