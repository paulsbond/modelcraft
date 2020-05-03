from typing import Iterator, Iterable, List, Optional, Union
from dataclasses import dataclass
import itertools
import gemmi
import numpy
import pandas


@dataclass
class ColumnRef:
    label: str
    dataset: str = ""
    crystal: str = ""
    project: str = ""

    def matching_columns(self, mtz: gemmi.Mtz) -> Iterator[gemmi.Mtz.Column]:
        for column in mtz.columns:
            if (
                self.label == column.label
                and self.dataset in ("", column.dataset.dataset_name)
                and self.crystal in ("", column.dataset.crystal_name)
                and self.project in ("", column.dataset.project_name)
            ):
                yield column

    def find_column(self, mtz: gemmi.Mtz) -> gemmi.Mtz.Column:
        matching = list(self.matching_columns(mtz))
        if len(matching) == 0:
            raise ValueError(f"No columns matching '{self}' in MTZ")
        if len(matching) > 1:
            raise ValueError(f"Multiple columns matching '{self}' in MTZ")
        return matching[0]

    def __str__(self):
        parts = (self.project, self.crystal, self.dataset, self.label)
        return "/".join(part for part in parts if part != "")


def expand_label(label: str) -> str:
    i = label.rindex(".") if "." in label else -1
    prefix = label[: i + 1]
    suffix = label[i + 1 :]
    if suffix[-4:] == "ABCD":
        return ",".join(f"{prefix}{suffix}.{x}" for x in ("A", "B", "C", "D"))
    if suffix[:2] == "HL" or suffix[-2:] == "HL":
        return ",".join(f"{prefix}{suffix}{x}" for x in ("A", "B", "C", "D"))
    if suffix in ("F_sigF", "I_sigI", "F_phi", "phi_fom"):
        return ",".join(f"{prefix}{suffix}.{x}" for x in suffix.split("_"))
    return label


def column_refs(columns: str) -> List[ColumnRef]:
    columns = columns.replace("*", "")
    split = [""] * 3 + columns.split("/")
    project, crystal, dataset, label = split[-4:]
    if "," not in label:
        label = expand_label(label)
    return [ColumnRef(label, dataset, crystal, project) for label in label.split(",")]


class DataItem(gemmi.Mtz):
    def __init__(self, mtz: gemmi.Mtz, columns: Union[str, Iterable[gemmi.Mtz.Column]]):
        super().__init__()
        self.cell = mtz.cell
        self.spacegroup = mtz.spacegroup
        self.add_dataset("HKL_base")
        if isinstance(columns, str):
            columns = [ref.find_column(mtz) for ref in column_refs(columns)]
        columns = list(mtz.columns[:3]) + list(columns)
        for column in columns:
            self.add_column(column.label, column.type)
        data = numpy.stack(columns, axis=1)
        self.set_data(data)
        # TODO: Add self.update_reso() when it is available

    def label(self, index: Optional[int] = None) -> str:
        if index is None:
            return ",".join(column.label for column in self.columns[3:])
        return self.columns[index + 3].label

    def data_frame(self) -> pandas.DataFrame:
        data = numpy.array(self, copy=False)
        return pandas.DataFrame(data=data, columns=self.column_labels())

    @classmethod
    def search(cls, mtz: gemmi.Mtz, sequential: bool = True):
        if cls is DataItem:
            raise TypeError("Cannot search for abstract data items")
        if sequential:
            mtz_types = [col.type for col in mtz.columns]
            num_cols = len(cls._types)
            i = 0
            while i < len(mtz_types) - num_cols + 1:
                if mtz_types[i : i + num_cols] == cls._types:
                    columns = mtz.columns[i : i + num_cols]
                    yield cls(mtz, columns)
                    i += num_cols
                else:
                    i += 1
        else:
            columns = []
            for column_type in cls._types:
                columns.append([col for col in mtz.columns if col.type == column_type])
            for combination in itertools.product(*columns):
                yield cls(mtz, combination)


class FsigF(DataItem):
    _types = ["F", "Q"]


class FreeRFlag(DataItem):
    _types = ["I"]


class FPhi(DataItem):
    _types = ["F", "P"]


class ABCD(DataItem):
    _types = ["A", "A", "A", "A"]


class PhiFom(DataItem):
    _types = ["P", "W"]


def _combine_data_items(items: List[DataItem]) -> gemmi.Mtz:
    column_labels = ["H", "K", "L"]
    column_types = ["H", "H", "H"]
    for item in items:
        for column in item.columns[3:]:
            column_labels.append(column.label)
            column_types.append(column.type)
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
    items = [item for item in items if item is not None]
    mtz = _combine_data_items(items)
    mtz.write_to_file(path)
