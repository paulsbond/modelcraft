from PySide2.QtWidgets import QApplication, QDialog
import os
import sys


class Window(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("ModelCraft")


def main():
    if "CCP4" not in os.environ:
        sys.exit("CCP4 environment not set")
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
