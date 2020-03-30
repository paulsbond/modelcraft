import gemmi
import itertools
import os


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

    def _data_items(self, mtz, types, sequential=True):
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
