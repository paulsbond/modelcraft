import os
import sys
from modelcraft.reflections import DataFile
from PySide2.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
)


class ReflectionGroupBox(QGroupBox):
    def __init__(self, parentDialog):
        super().__init__("Reflection Data")
        self.parentDialog = parentDialog

        layout = QVBoxLayout()
        self.setLayout(layout)

        hbox = QHBoxLayout()
        button = QPushButton("Select File")
        button.clicked.connect(self.select_file)
        hbox.addWidget(button)
        self.file_label = QLabel("No file selected.")
        hbox.addWidget(self.file_label)
        layout.addLayout(hbox)

        self.spacegroup_label = QLabel("Spacegroup: ")
        layout.addWidget(self.spacegroup_label)

    def select_file(self):
        path, _ = QFileDialog.getOpenFileName(
            parent=self.parentDialog,
            caption="Select Reflection Data File",
            filter="MTZ Files [*.mtz] (*.mtz) ;; All Files [*.*] (*.*)",
            options=QFileDialog.Options(),
        )
        if path != "":
            self.file_label.setText(os.path.basename(path))
            hklin = DataFile(path)
            self.spacegroup_label.setText("Spacegroup: " + hklin.spacegroup.hm)


class Window(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("ModelCraft")
        layout = QVBoxLayout()
        self.setLayout(layout)
        layout.addWidget(ReflectionGroupBox(self))


def main():
    if "CCP4" not in os.environ:
        sys.exit("CCP4 environment not set")
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
