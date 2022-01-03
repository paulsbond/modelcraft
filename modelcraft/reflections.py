from typing import Iterator, Iterable, List, Optional, Union
import itertools
import re
import gemmi
import numpy
import pandas


class ColumnRef:
    def __init__(
        self, label: str, dataset: str = "", crystal: str = "", project: str = ""
    ):
        self.label = label
        self.dataset = dataset
        self.crystal = crystal
        self.project = project

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


def contract_label(label: str) -> str:
    for name, pattern in (
        ("ABCD", r"(.*)ABCD\.A,(.*)ABCD\.B,(.*)ABCD\.C,(.*)ABCD\.D"),
        ("F_phi", r"(.*)F_phi\.F,(.*)F_phi\.phi"),
        ("F_sigF", r"(.*)F_sigF\.F,(.*)F_sigF\.sigF"),
        ("I_sigI", r"(.*)I_sigI\.I,(.*)I_sigI\.sigI"),
        ("phi_fom", r"(.*)phi_fom\.phi,(.*)phi_fom\.fom"),
    ):
        match = re.match(pattern, label)
        if match and len(set(match.groups())) == 1:
            return match.group(1) + name
    return label


def column_refs(columns: str) -> List[ColumnRef]:
    columns = columns.replace("*", "")
    columns = columns.rstrip("/")
    split = [""] * 3 + columns.split("/")
    project, crystal, dataset, label = split[-4:]
    label = label.lstrip("[")
    label = label.rstrip("]")
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
        self.types = "".join(col.type for col in columns)
        columns = list(mtz.columns[:3]) + list(columns)
        for column in columns:
            self.add_column(column.label, column.type)
        data = numpy.stack(columns, axis=1)
        data = data[~numpy.isnan(data[:, 3:]).all(axis=1)]
        self.set_data(data)
        self.update_reso()

    def label(self, index: Optional[int] = None) -> str:
        if index is None:
            return ",".join(column.label for column in self.columns[3:])
        return self.columns[index + 3].label

    def data_frame(self, copy=False) -> pandas.DataFrame:
        data = numpy.array(self, copy=copy)
        return pandas.DataFrame(data=data, columns=self.column_labels())

    @classmethod
    def search(cls, mtz: gemmi.Mtz, types: str, sequential: bool = True):
        types = list(types)
        if sequential:
            mtz_types = [col.type for col in mtz.columns]
            num_cols = len(types)
            i = 0
            while i < len(mtz_types) - num_cols + 1:
                if mtz_types[i : i + num_cols] == types:
                    columns = mtz.columns[i : i + num_cols]
                    yield cls(mtz, columns)
                    i += num_cols
                else:
                    i += 1
        else:
            columns = []
            for column_type in types:
                columns.append([col for col in mtz.columns if col.type == column_type])
            for combination in itertools.product(*columns):
                yield cls(mtz, combination)


def write_mtz(
    path: str, items: List[DataItem], labels: Optional[List[str]] = None
) -> None:
    labels = labels or [None] * len(items)
    column_labels = ["H", "K", "L"]
    column_types = ["H", "H", "H"]
    data = None
    for item, label in zip(items, labels):
        if item is not None:
            if data is None:
                data = item.data_frame()
            else:
                data = data.merge(item.data_frame(), on=["H", "K", "L"], how="outer")
            column_labels.extend(
                (col.label for col in item.columns[3:])
                if label is None
                else label.split(",")
            )
            column_types.extend(col.type for col in item.columns[3:])
    mtz = gemmi.Mtz()
    mtz.cell = items[0].cell
    mtz.spacegroup = items[0].spacegroup
    mtz.add_dataset("HKL_base")
    for label, type_ in zip(column_labels, column_types):
        mtz.add_column(label, type_)
    mtz.set_data(data.to_numpy())
    mtz.write_to_file(path)
