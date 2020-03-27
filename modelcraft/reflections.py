import os
import gemmi


class ReflectionFile:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        mtz = gemmi.read_mtz_file(path)

        self.cell = mtz.cell
        self.spacegroup = mtz.spacegroup
        self.num_reflections = mtz.nreflections
        self.resolution_high = mtz.resolution_high()
        self.resolution_low = mtz.resolution_low()

        self.columns = {col.label for col in mtz.columns}
        self.fsigfs = columns_with_types(mtz, ["F", "Q"])
        self.frees = columns_with_types(mtz, ["I"])
        self.abcds = columns_with_types(mtz, ["A", "A", "A", "A"])
        self.phifoms = set()
        self.fphis = set()
        for i in range(len(mtz.columns)):
            col1 = mtz.columns[0]
            if col1.type == "F" and i < len(mtz.columns) - 1:
                col2 = mtz.columns[i + 1]
                if col2.type == "Q":
                    labels = ",".join((col1.label, col2.label))
                    self.fsigfs.add(labels)
                    i + 1
            elif col1.type == "I":
                self.frees.add(col1.label)

    def __contains__(self, column):
        return column in self.columns


def columns_with_types(mtz, types):
    pass


def fo_columns(mtz):
    def differ_by_sig(label1, label2):
        label1 = label1.lower().replace("sig", "")
        label2 = label2.lower().replace("sig", "")
        return label1 == label2

    for f in mtz.columns_with_type("F"):
        for sig in mtz.columns_with_type("Q"):
            if differ_by_sig(f.label, sig.label):
                yield "%s,%s" % (f.label, sig.label)


def free_columns(mtz):
    for column in mtz.columns_with_type("I"):
        if "free" in column.label.lower():
            yield column.label


def hl_columns(mtz):
    def belong_to_the_same_group(label1, label2):
        if len(label1) != len(label2):
            return False
        differences = 0
        for i in range(len(label1)):
            if label1[i] != label2[i]:
                differences += 1
        return differences == 1

    groups = []
    for column in mtz.columns_with_type("A"):
        grouped = False
        for group in groups:
            if belong_to_the_same_group(column.label, group[0]):
                grouped = True
                group.append(column.label)
                if len(group) == 4:
                    yield ",".join(sorted(group))
        if not grouped:
            groups.append([column.label])


def phifom_columns(mtz):
    def is_phi(column):
        if column.type != "P":
            return False
        label = column.label.lower()
        return "ph" in label and "del" not in label

    def is_fom(column):
        if column.type != "W":
            return False
        label = column.label.lower()
        return "fom" in label

    phis = [column for column in mtz.columns if is_phi(column)]
    foms = [column for column in mtz.columns if is_fom(column)]
    for phi in phis:
        for fom in foms:
            yield "%s,%s" % (phi.label, fom.label)
