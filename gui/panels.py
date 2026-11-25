"""Display widgets: the verdict card, the factor breakdown, and the structure view."""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QFont, QPixmap
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from gui.depict import LEGEND, render_png, rgb_css
from pipeline import ScanResult

# Severity colour per potency category. Category 1 is the most potent.
CATEGORY_COLOUR = {
    1: "#B3261E",
    2: "#C2410C",
    3: "#B45309",
    4: "#3F6212",
    5: "#3F6212",
}
NEUTRAL = "#44403C"


def pretty(text: str) -> str:
    """The core keeps factor names in ASCII so console output is safe on Windows;
    the window can afford the real Greek letters."""
    return text.replace("alpha-", "α-").replace("beta-", "β-")


class VerdictCard(QFrame):
    """The headline answer: category, acceptable intake, and where it came from."""

    def __init__(self) -> None:
        super().__init__()
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setMinimumHeight(132)

        self._category = QLabel("")
        self._category.setFont(QFont("Segoe UI", 20, QFont.Weight.Bold))

        self._limit = QLabel("Enter a structure, name or CAS number to begin.")
        self._limit.setFont(QFont("Segoe UI", 15))
        self._limit.setWordWrap(True)

        self._basis = QLabel("")
        self._basis.setWordWrap(True)
        self._basis.setFont(QFont("Segoe UI", 9))

        self._ppm = QLabel("")
        self._ppm.setFont(QFont("Segoe UI", 10))

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 12, 16, 12)
        layout.setSpacing(3)
        for w in (self._category, self._limit, self._basis, self._ppm):
            layout.addWidget(w)
        layout.addStretch(1)
        self._paint(NEUTRAL)

    def _paint(self, colour: str) -> None:
        self.setStyleSheet(
            f"QFrame {{ background: {colour}; border-radius: 8px; }}"
            "QLabel { color: white; background: transparent; }")

    def clear(self) -> None:
        self._category.setText("")
        self._limit.setText("Working...")
        self._basis.setText("")
        self._ppm.setText("")
        self._paint(NEUTRAL)

    def show_error(self, message: str) -> None:
        self._category.setText("Not assessed")
        self._limit.setText(message)
        self._basis.setText("")
        self._ppm.setText("")
        self._paint(NEUTRAL)

    def show(self, result: ScanResult, mdd_mg: float | None = None) -> None:
        if result.error:
            self.show_error(result.error)
            return

        ai = result.final_ai_ng_day
        cat = result.final_category

        if ai is None:
            reason = (result.cpca.not_applicable_reason
                      if result.cpca else "No limit could be determined.")
            self.show_error(reason or "CPCA does not apply to this structure.")
            return

        self._category.setText(f"Potency Category {cat}" if cat else "Limit established")
        self._limit.setText(f"Acceptable intake  {ai:g} ng/day")
        self._paint(CATEGORY_COLOUR.get(cat or 0, NEUTRAL))

        basis = result.final_basis
        if result.overridden and result.cpca:
            basis += (f"   ·   overrides the CPCA prediction of Category "
                      f"{result.cpca.category} ({result.cpca.ai_ng_day:g} ng/day)")
        if result.published and result.published.is_interim:
            basis += "   ·   INTERIM, confirm it is still current"
        self._basis.setText(basis)

        if mdd_mg:
            ppm = result.control_limit_ppm(mdd_mg)
            self._ppm.setText(
                f"At a max daily dose of {mdd_mg:g} mg  →  control limit "
                f"{ppm:.4g} ppm" if ppm else "")
        else:
            self._ppm.setText("")


class FactorTable(QTableWidget):
    """Every line item that contributed to the potency score."""

    HEADERS = ["Score", "Factor", "Source", ""]

    def __init__(self) -> None:
        super().__init__(0, len(self.HEADERS))
        self.setHorizontalHeaderLabels(self.HEADERS)
        self.verticalHeader().setVisible(False)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setAlternatingRowColors(True)
        self.setWordWrap(True)

        header = self.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)

    def show_result(self, result: ScanResult) -> None:
        self.setRowCount(0)
        cpca = result.cpca
        if cpca is None or not cpca.applicable:
            return
        site = cpca.governing_site
        if site is None:
            return

        rows: list[tuple[str, str, str, str]] = []

        if site.alpha_h_counts:
            a, b = site.alpha_h_counts
            rows.append(("", f"α-hydrogen counts (lowest first): {a},{b}",
                         "Figure 1", ""))

        if site.triage_exit:
            rows.append(("→ 5", pretty(site.triage_exit), "Figure 1", ""))
        else:
            for f in site.factors:
                label = pretty(f.name) + (f":  {pretty(f.detail)}" if f.detail else "")
                rows.append((f"{f.score:+d}", label, f.reference,
                             "review" if f.uncertain else ""))
            rows.append((f"= {site.potency_score}", "Potency Score", "Appendix A", ""))

        for note in site.notes:
            rows.append(("", pretty(note), "", "review"))

        self.setRowCount(len(rows))
        for r, (score, name, ref, flag) in enumerate(rows):
            score_item = QTableWidgetItem(score)
            score_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            if score.startswith("+"):
                score_item.setForeground(QColor("#15803D"))
            elif score.startswith("-"):
                score_item.setForeground(QColor("#B3261E"))
            if score.startswith("="):
                font = score_item.font()
                font.setBold(True)
                score_item.setFont(font)

            name_item = QTableWidgetItem(name)
            if score.startswith("="):
                font = name_item.font()
                font.setBold(True)
                name_item.setFont(font)

            flag_item = QTableWidgetItem("⚠ review" if flag else "")
            if flag:
                flag_item.setForeground(QColor("#B45309"))
                flag_item.setToolTip(
                    "FDA defers to Cross & Ponting (2021) rather than listing these "
                    "groups. Have a qualified assessor confirm this match.")

            for col, item in enumerate((score_item, name_item,
                                        QTableWidgetItem(ref), flag_item)):
                self.setItem(r, col, item)

        self.resizeRowsToContents()


PANEL_STYLE = ("QLabel { background: white; border: 1px solid #D6D3D1; "
               "border-radius: 6px; color: #78716C; }")
CAPTION_STYLE = "color: #57534E; font-size: 11px; border: none; background: transparent;"


class StructureView(QWidget):
    """The molecule, with the scored atoms colour-coded.

    When the input was a picture, the source picture is shown beside our redrawing of
    what was read from it. Recognition gets roughly one structure in five wrong, so the
    only safe design is to put both in front of the user and let them compare.
    """

    def __init__(self) -> None:
        super().__init__()
        self._source_pixmap: QPixmap | None = None
        self._drawn_pixmap: QPixmap | None = None

        self._source_caption = QLabel("The picture you gave me")
        self._source_caption.setStyleSheet(CAPTION_STYLE)
        self._source = QLabel()
        self._source.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._source.setMinimumSize(210, 260)
        self._source.setScaledContents(False)
        self._source.setStyleSheet(PANEL_STYLE)

        self._source_box = QWidget()
        source_layout = QVBoxLayout(self._source_box)
        source_layout.setContentsMargins(0, 0, 0, 0)
        source_layout.setSpacing(3)
        source_layout.addWidget(self._source_caption)
        source_layout.addWidget(self._source, 1)
        self._source_box.setVisible(False)

        self._drawn_caption = QLabel("")
        self._drawn_caption.setStyleSheet(CAPTION_STYLE)
        self._drawn_caption.setVisible(False)
        self._image = QLabel("Drop a structure image here, or enter a SMILES string")
        self._image.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._image.setWordWrap(True)
        self._image.setMinimumSize(210, 260)
        self._image.setStyleSheet(PANEL_STYLE)

        drawn_box = QWidget()
        drawn_layout = QVBoxLayout(drawn_box)
        drawn_layout.setContentsMargins(0, 0, 0, 0)
        drawn_layout.setSpacing(3)
        drawn_layout.addWidget(self._drawn_caption)
        drawn_layout.addWidget(self._image, 1)

        panels = QHBoxLayout()
        panels.setSpacing(10)
        panels.addWidget(self._source_box, 1)
        panels.addWidget(drawn_box, 1)

        self._warning = QLabel(
            "Check the two drawings match. Recognition reads roughly four structures "
            "in five correctly, so the one on the right is a reading, not a fact.")
        self._warning.setWordWrap(True)
        self._warning.setStyleSheet(
            "QLabel { background: #FEF3C7; border: 1px solid #F59E0B; "
            "border-radius: 4px; padding: 5px; color: #78350F; font-size: 11px; }")
        self._warning.setVisible(False)

        legend = QHBoxLayout()
        legend.setSpacing(12)
        for label, colour in LEGEND:
            swatch = QLabel("●")
            swatch.setStyleSheet(f"color: {rgb_css(colour)}; font-size: 15px;")
            text = QLabel(label)
            text.setStyleSheet(CAPTION_STYLE)
            legend.addWidget(swatch)
            legend.addWidget(text)
        legend.addStretch(1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(panels, 1)
        layout.addWidget(self._warning)
        layout.addLayout(legend)

    # state

    def clear(self) -> None:
        self._source_pixmap = None
        self._drawn_pixmap = None
        self._image.setPixmap(QPixmap())
        self._image.setText("Drop a structure image here, or enter a SMILES string")
        self._drawn_caption.setVisible(False)
        self._source_box.setVisible(False)
        self._warning.setVisible(False)

    @staticmethod
    def _fit(label: QLabel, pixmap: QPixmap | None) -> None:
        """QLabel clips a pixmap larger than itself rather than scaling it, so the
        original is kept and refitted whenever the panel changes size."""
        if pixmap is None or pixmap.isNull():
            return
        width = max(label.width() - 8, 120)
        height = max(label.height() - 8, 120)
        label.setPixmap(pixmap.scaled(
            width, height,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation))
        label.setText("")

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._fit(self._source, self._source_pixmap)
        self._fit(self._image, self._drawn_pixmap)

    def show_source_image(self, path) -> None:
        """Display the picture on the left while recognition runs."""
        pixmap = QPixmap(str(path))
        if pixmap.isNull():
            return
        self._source_pixmap = pixmap
        self._source_box.setVisible(True)
        self._drawn_caption.setText("What it read")
        self._drawn_caption.setVisible(True)
        self._warning.setVisible(True)
        self._drawn_pixmap = None
        self._image.setPixmap(QPixmap())
        self._image.setText("Reading the structure...")
        self._fit(self._source, pixmap)

    def show_smiles(self, smiles: str | None) -> None:
        if not smiles:
            self._drawn_pixmap = None
            self._image.setPixmap(QPixmap())
            self._image.setText("No structure")
            return
        png = render_png(smiles, 620, 450)
        if png is None:
            self._drawn_pixmap = None
            self._image.setPixmap(QPixmap())
            self._image.setText("Could not draw this structure")
            return
        pixmap = QPixmap()
        pixmap.loadFromData(png)
        self._drawn_pixmap = pixmap
        self._fit(self._image, pixmap)
