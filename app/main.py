import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QComboBox, QStackedWidget, QScrollArea, QSizePolicy, QPushButton,
    QButtonGroup, QFrame
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor, QPalette

from brdf_form import BRDFForm
from mie_form import MieForm
from reflect_form import ReflectForm
from rcw_form import RCWForm


AFFILIATION = "Woojun Lee · Purdue University · School of Mechanical Engineering"


def apply_light_theme(app):

    app.setStyle("Fusion")
    palette = QPalette()
    colors = {
        QPalette.Window: "#ffffff",
        QPalette.WindowText: "#202124",
        QPalette.Base: "#ffffff",
        QPalette.AlternateBase: "#f5f5f5",
        QPalette.ToolTipBase: "#ffffff",
        QPalette.ToolTipText: "#202124",
        QPalette.Text: "#202124",
        QPalette.Button: "#f2f2f2",
        QPalette.ButtonText: "#202124",
        QPalette.BrightText: "#ffffff",
        QPalette.Light: "#ffffff",
        QPalette.Midlight: "#f5f5f5",
        QPalette.Mid: "#c7c7c7",
        QPalette.Dark: "#999999",
        QPalette.Shadow: "#666666",
        QPalette.Highlight: "#2673c9",
        QPalette.HighlightedText: "#ffffff",
        QPalette.Link: "#1565c0",
        QPalette.LinkVisited: "#6a3dad",
        QPalette.PlaceholderText: "#757575",
    }
    for role, color in colors.items():
        palette.setColor(role, QColor(color))
    for role in (QPalette.WindowText, QPalette.Text, QPalette.ButtonText):
        palette.setColor(QPalette.Disabled, role, QColor("#808080"))
    palette.setColor(QPalette.Disabled, QPalette.Highlight, QColor("#e0e0e0"))
    palette.setColor(QPalette.Disabled, QPalette.HighlightedText, QColor("#808080"))
    app.setPalette(palette)


class SCATMECHGui(QWidget):
    def closeEvent(self, event):
        for form in (self.brdf_form, self.rcw_form):
            viewer = getattr(form, "_geometry_viewer", None)
            if viewer is not None:
                viewer.close()
        super().closeEvent(event)

    def __init__(self):
        apply_light_theme(QApplication.instance())
        super().__init__()
        self.setWindowTitle("SCATMECH GUI")


        root = QVBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)


        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.setSpacing(6)

        self.btn_group = QButtonGroup(self)
        self.btn_group.setExclusive(True)

        def make_btn(text):
            b = QPushButton(text, self)
            b.setCheckable(True)
            b.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            return b

        self.btn_brdf    = make_btn("BRDF")
        self.btn_rcw     = make_btn("RCW")
        self.btn_reflect = make_btn("Reflect")
        self.btn_mie     = make_btn("Mie")

        buttons = [self.btn_brdf, self.btn_rcw, self.btn_reflect, self.btn_mie]
        for i, b in enumerate(buttons):
            self.btn_group.addButton(b, i)
            top.addWidget(b)
            top.setStretch(i, 1)

        root.addLayout(top)


        self.stack = QStackedWidget(self)


        self.brdf_form    = BRDFForm()
        self.reflect_form = ReflectForm()
        self.mie_form     = MieForm()
        self.rcw_form = RCWForm()

        def wrap_in_scroll(w):

            form_widget = w.main_layout.takeAt(0).widget()
            sa = QScrollArea(w)
            sa.setWidgetResizable(True)
            sa.setWidget(form_widget)
            sa.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Expanding)
            sa.setMinimumWidth(300)
            w.main_layout.insertWidget(0, sa, 1)
            w.main_layout.setStretch(1, 1)
            return w


        self.stack.addWidget(wrap_in_scroll(self.brdf_form))
        self.stack.addWidget(wrap_in_scroll(self.rcw_form))
        self.stack.addWidget(wrap_in_scroll(self.reflect_form))
        self.stack.addWidget(wrap_in_scroll(self.mie_form))

        root.addWidget(self.stack, 1)


        footer = QFrame(self)
        footer.setObjectName("footer")
        footer.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        footer.setMinimumHeight(24)
        footer.setStyleSheet(
            "#footer { border-top: 1px solid #e5e5e5; }"
            "#footer QLabel { color: #666; font-size: 11px; }"
        )
        footer_layout = QHBoxLayout(footer)
        footer_layout.setContentsMargins(8, 2, 8, 6)
        footer_layout.setSpacing(0)

        self.affiliation_label = QLabel(AFFILIATION, footer)
        self.affiliation_label.setAlignment(Qt.AlignCenter)

        footer_layout.addStretch(1)
        footer_layout.addWidget(self.affiliation_label, 0, Qt.AlignCenter)
        footer_layout.addStretch(1)

        root.addWidget(footer, 0)


        def set_page(idx):
            self.stack.setCurrentIndex(idx)
            btn = self.btn_group.button(idx)
            if btn and not btn.isChecked():
                btn.setChecked(True)

        self.btn_brdf.clicked.connect(lambda: set_page(0))
        self.btn_rcw.clicked.connect(lambda: set_page(1))
        self.btn_reflect.clicked.connect(lambda: set_page(2))
        self.btn_mie.clicked.connect(lambda: set_page(3))


        set_page(0)
        self.btn_brdf.setChecked(True)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    gui = SCATMECHGui()
    gui.show()
    sys.exit(app.exec_())
