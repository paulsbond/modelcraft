import itertools
import gemmi
import numpy
import pandas


def matching_columns(mtz, pattern):
    pattern = pattern.replace("*", "")
    split = [""] * 3 + pattern.split("/")
    project, crystal, dataset, label = split[-4:]
    for column in mtz.columns:
        if (
            label in ("", column.label)
            and dataset in ("", column.dataset.dataset_name)
            and crystal in ("", column.dataset.crystal_name)
            and project in ("", column.dataset.project_name)
        ):
            yield column


def find_column(mtz, pattern):
    matching = list(matching_columns(mtz, pattern))
    if len(matching) == 0:
        raise ValueError(f"No columns matching '{pattern}' in MTZ")
    if len(matching) > 1:
        raise ValueError(f"Multiple columns matching '{pattern}' in MTZ")
    return matching[0]


def find_columns(mtz, pattern):
    split = pattern.split(",")
    return [find_column(mtz, pattern) for pattern in split]


class DataItem(gemmi.Mtz):
    def __init__(self, mtz, columns):
        super().__init__()
        if isinstance(columns, str):
            columns = find_columns(mtz, columns)
        self.cell = mtz.cell
        self.spacegroup = mtz.spacegroup
        self.add_dataset("HKL_base")
        columns = list(mtz.columns[:3]) + list(columns)
        for column in columns:
            self.add_column(column.label, column.type)
        data = numpy.stack(columns, axis=1)
        self.set_data(data)

    def label(self, index=None):
        if index is None:
            return ",".join(column.label for column in self.columns[3:])
        return self.columns[index + 3].label

    def data_frame(self):
        data = numpy.array(self, copy=False)
        return pandas.DataFrame(data=data, columns=self.column_labels())

    @classmethod
    def search(cls, mtz, sequential=True):
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


def _combine_data_items(items):
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


def write_mtz(path, items):
    mtz = _combine_data_items(items)
    mtz.write_to_file(path)
