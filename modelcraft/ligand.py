from typing import Optional


class Ligand:
    def __init__(self, code: str, copies: Optional[int] = None):
        self.code = code
        self.copies = copies

    def __eq__(self, other) -> bool:
        if isinstance(other, Ligand):
            return self.code == other.code
        return NotImplemented

    @classmethod
    def from_component_json(cls, component: dict) -> "Ligand":
        return cls(code=component["code"], copies=component.get("copies"))

    @classmethod
    def from_pdbe_molecule_dict(cls, mol: dict) -> "Ligand":
        return cls(code=mol["chem_comp_ids"][0], copies=mol["number_of_copies"])

    def to_component_json(self) -> dict:
        return {"code": self.code, "copies": self.copies}
