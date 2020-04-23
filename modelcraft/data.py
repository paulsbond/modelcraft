from typing import List
import itertools
import gemmi
import numpy
import pandas
from modelcraft.utils import all_equal, contains_duplicates


def cells_are_equal(*cells: gemmi.UnitCell):
    return (
        all_equal(cell.a for cell in cells)
        and all_equal(cell.b for cell in cells)
        and all_equal(cell.c for cell in cells)
        and all_equal(cell.alpha for cell in cells)
        and all_equal(cell.beta for cell in cells)
        and all_equal(cell.gamma for cell in cells)
    )


def find_column(mtz: gemmi.Mtz, label: str):
    # TODO: Support project/crystal/dataset/column names
    if label not in mtz.column_labels():
        raise ValueError(f"No '{label}' column in MTZ")
    if mtz.column_labels().count(label) > 1:
        raise ValueError(f"Multiple '{label}' columns in MTZ")
    return mtz.column_with_label(label)


class DataItem(gemmi.Mtz):
    def __init__(self, mtz: gemmi.Mtz, labels: List[str]):
        super().__init__()
        self.cell = mtz.cell
        self.spacegroup = mtz.spacegroup
        self.add_dataset("HKL_base")
        columns = []
        for label in ["H", "K", "L"] + labels:
            column = find_column(mtz, label)
            if not cells_are_equal(column.dataset.cell, mtz.cell):
                raise ValueError("Dataset cell not equal to the MTZ cell")
            columns.append(column)
            self.add_column(column.label, column.type)
        data = numpy.stack(columns, axis=1)
        self.set_data(data)

    def label(self, index: int = None) -> str:
        if index is None:
            return ",".join(column.label for column in self.columns[3:])
        return self.columns[index + 3].label

    def data_frame(self) -> pandas.DataFrame:
        data = numpy.array(self, copy=False)
        return pandas.DataFrame(data=data, columns=self.column_labels())


def _combine_data_items(items: List[DataItem]) -> gemmi.Mtz:
    assert len(items) > 1
    assert cells_are_equal(*(item.cell for item in items))
    assert all_equal(item.spacegroup for item in items)
    column_labels = ["H", "K", "L"]
    column_types = ["H", "H", "H"]
    for item in items:
        for column in item.columns[3:]:
            column_labels.append(column.label)
            column_types.append(column.type)
    assert not contains_duplicates(column_labels)
    data_frame = items[0].data_frame()
    for item in items[1:]:
        data_frame = pandas.merge(data_frame, item.data_frame(), on=["H", "K", "L"])
    mtz = gemmi.Mtz()
    mtz.cell = items[0].cell
    mtz.spacegroup = items[0].spacegroup
    mtz.add_dataset("HKL_base")
    for i, label in enumerate(column_labels):
        mtz.add_column(label, column_types[i])
    mtz.set_data(data_frame.to_numpy())
    return mtz


def write_mtz(path: str, items: List[DataItem]) -> None:
    mtz = _combine_data_items(items)
    mtz.write_to_file(path)


# TODO: Return labels used to uniquely identify the columns
def _columns(
    mtz: gemmi.Mtz, types: List[str], sequential: bool = True
) -> List[gemmi.Mtz.Column]:
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
