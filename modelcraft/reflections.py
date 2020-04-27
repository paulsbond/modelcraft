from abc import ABC
from typing import List, Type, Tuple
from dataclasses import dataclass
import itertools
import gemmi
import numpy
import pandas
from modelcraft.utils import all_equal, contains_duplicates


def find_column(mtz: gemmi.Mtz, pattern: str) -> gemmi.Mtz.Column:
    pattern = pattern.replace("*", "")
    split = [""] * 3 + pattern.split("/")
    project, crystal, dataset, label = split[-4:]
    matching = []
    for column in mtz.columns:
        if (
            label in ("", column.label)
            and dataset in ("", column.dataset.dataset_name)
            and crystal in ("", column.dataset.crystal_name)
            and project in ("", column.dataset.project_name)
        ):
            matching.append(column)
    if len(matching) == 0:
        raise ValueError(f"No columns matching '{pattern}' in MTZ")
    if len(matching) > 1:
        raise ValueError(f"Multiple columns matching '{pattern}' in MTZ")
    return matching[0]


# class DataItem(ABC):
#     def __init__(self, mtz: gemmi.Mtz, columns: List[str]):
#         self._mtz = gemmi.Mtz()
#         self._mtz.cell = mtz.cell
#         self._mtz.spacegroup = mtz.spacegroup
#         self._mtz.add_dataset("HKL_base")

#         for label in ["H", "K", "L"] + labels:
#             column = find_column(mtz, label)
#             if not cells_are_equal(column.dataset.cell, mtz.cell):
#                 raise ValueError("Dataset cell not equal to the MTZ cell")
#             columns.append(column)
#             self._mtz.add_column(column.label, column.type)
#         data = numpy.stack(columns, axis=1)
#         self._mtz.set_data(data)

#     def label(self, index: int = None) -> str:
#         if index is None:
#             return ",".join(column.label for column in self._mtz.columns[3:])
#         return self._mtz.columns[index + 3].label

#     def data_frame(self) -> pandas.DataFrame:
#         data = numpy.array(self, copy=False)
#         return pandas.DataFrame(data=data, columns=self._mtz.column_labels())

#     @classmethod
#     def from_labels(cls, mtz: gemmi.Mtz, labels: List[str]) -> cls:
#         columns = [find_column(mtz, label) for label in labels]


# class FsigF(DataItem):
#     _types = ["F", "Q"]
#     _sequential = True


# def _combine_data_items(items: List[DataItem]) -> gemmi.Mtz:
#     assert len(items) > 1
#     assert cells_are_equal(*(item.cell for item in items))
#     assert all_equal(item.spacegroup for item in items)
#     column_labels = ["H", "K", "L"]
#     column_types = ["H", "H", "H"]
#     for item in items:
#         for column in item.columns[3:]:
#             column_labels.append(column.label)
#             column_types.append(column.type)
#     assert not contains_duplicates(column_labels)
#     data_frame = items[0].data_frame()
#     for item in items[1:]:
#         data_frame = pandas.merge(data_frame, item.data_frame(), on=["H", "K", "L"])
#     mtz = gemmi.Mtz()
#     mtz.cell = items[0].cell
#     mtz.spacegroup = items[0].spacegroup
#     mtz.add_dataset("HKL_base")
#     for i, label in enumerate(column_labels):
#         mtz.add_column(label, column_types[i])
#     mtz.set_data(data_frame.to_numpy())
#     return mtz


# def write_mtz(path: str, items: List[DataItem]) -> None:
#     mtz = _combine_data_items(items)
#     mtz.write_to_file(path)


# # TODO: Return labels used to uniquely identify the columns
# def _columns(
#     mtz: gemmi.Mtz, types: List[str], sequential: bool = True
# ) -> List[gemmi.Mtz.Column]:
#     items = set()
#     if sequential:
#         mtz_types = [col.type for col in mtz.columns]
#         i = 0
#         while i < len(mtz_types) - len(types) + 1:
#             if mtz_types[i : i + len(types)] == types:
#                 item = ",".join(col.label for col in mtz.columns[i : i + len(types)])
#                 items.add(item)
#                 i += 3
#             i += 1
#     else:
#         columns = [mtz.columns_with_type(t) for t in types]
#         for combination in itertools.product(*columns):
#             items.add(",".join(col.label for col in combination))
#     return items
