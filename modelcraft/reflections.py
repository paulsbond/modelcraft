import os
import gemmi


class ReflectionFile:
    def __init__(self, path, fsigf=None, free=None, abcd=None, phifom=None, fphi=None):
        self.path = os.path.abspath(path)
        self.fsigf = fsigf
        self.free = free
        self.abcd = abcd
        self.phifom = phifom
        self.fphi = fphi

    def resolution(self):
        mtz = gemmi.read_mtz_file(self.path)
        return round(mtz.resolution_high(), 2)


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
