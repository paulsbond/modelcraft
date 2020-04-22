from typing import List
import itertools
import os
import gemmi
import pandas


class MiniMtz(gemmi.Mtz):
    def __init__(self, mtz: gemmi.Mtz, *labels: str):
        super().__init__()
        self.cell = mtz.cell
        self.spacegroup = mtz.spacegroup
        self.label = ",".join(labels)
        self.columns = []
        for label in ("H", "K", "L") + labels:
            # TODO: throw error if dataset cell and spacegroup are not the same
            if label not in mtz.column_labels():
                raise ValueError(f"No '{label}' column in MTZ")
            if mtz.column_labels().count(label) > 1:
                raise ValueError(f"Multiple '{label}' columns in MTZ")
            self.columns.append(mtz.column_with_label(label))


def write_mtz(path: str, *data_items: DataItem) -> None:
    mtz = gemmi.Mtz()
    mtz.cell = data_items[0].mtz.cell
    mtz.spacegroup = data_items[0].mtz.spacegroup

    for data_item in data_items:
        for column in data_item.columns:
            mtz.add_column(column.label)
    pass
    # TODO: Throw an error if the column labels aren't unique
    # TODO: Create a combined MTZ file and write it out


def _columns(
    mtz: gemmi.Mtz, types: List[str], sequential: bool = True
) -> gemmi.MtzColumns:
    items = set()
    if sequential:
        mtz_types = [col.type for col in mtz.columns]
        i = 0
        while i < len(mtz_types) - len(types) + 1:
            if mtz_types[i : i + len(types)] == types:
                item = ",".join(col.label for col in mtz.columns[i : i + len(types)])
                items.add(item)
                i += 3
            i += 1
    else:
        columns = [mtz.columns_with_type(t) for t in types]
        for combination in itertools.product(*columns):
            items.add(",".join(col.label for col in combination))
    return items


class DataFile:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        mtz = gemmi.read_mtz_file(path)

        self.cell = mtz.cell
        self.spacegroup = mtz.spacegroup
        self.num_reflections = mtz.nreflections
        self.resolution_high = mtz.resolution_high()
        self.resolution_low = mtz.resolution_low()

        self.columns = {col.label for col in mtz.columns}
        self.fsigfs = self._data_items(mtz, ["F", "Q"])
        self.frees = self._data_items(mtz, ["I"])
        self.abcds = self._data_items(mtz, ["A", "A", "A", "A"])
        self.phifoms = self._data_items(mtz, ["P", "W"], sequential=False)
        self.fphis = self._data_items(mtz, ["F", "P"])

        self.fsigf = None
        self.free = None
        self.abcd = None
        self.phifom = None
        self.fwphiw = None
        self.fcphic = None

    def __contains__(self, columns):
        return all(col in self.columns for col in columns.split(","))
