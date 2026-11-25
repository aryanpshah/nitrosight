"""Desktop front end for the nitrosamine impurity scanner."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QGuiApplication
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

ROOT = Path(__file__).parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gui.panels import FactorTable, StructureView, VerdictCard  # noqa: E402
from gui.workers import BatchWorker, RecogniseWorker, ScanWorker  # noqa: E402
from osr import OsrService, is_available as osr_available  # noqa: E402
from pipeline import ScanResult, format_scan_report  # noqa: E402
from resolve import FDALimits, PubChemClient  # noqa: E402

KINDS = {"Auto-detect": "auto", "SMILES": "smiles", "Name": "name", "CAS number": "cas"}
IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".bmp", ".gif", ".tif", ".tiff"}

# FDA published a guidance document and a reference program that do not agree with each
# other on four rules. The default follows the document where it is clear and the
# program contradicts it, and follows the program where the document is ambiguous.
RULESETS = {
    "Balanced (default)": "balanced",
    "FDA guidance document": "guidance",
    "FDA reference program": "reference",
}

RULESET_HELP = (
    "FDA's guidance document and FDA's own reference program disagree on four rules.\n"
    "\n"
    "Balanced (default): follow the document where it is clear and the program\n"
    "contradicts it, follow the program where the document is ambiguous.\n"
    "\n"
    "FDA guidance document: follow the published document throughout.\n"
    "\n"
    "FDA reference program: reproduce FDA's featurize-nitrosamines tool throughout,\n"
    "including where it departs from the document.\n"
    "\n"
    "It changes the answer for about 1 in 80 structures. See\n"
    "docs/reference-implementation-divergences.md")

EXAMPLES = [
    ("NDEA", "CCN(N=O)CC"),
    ("N-nitroso-ketamine", "CN(N=O)C1(c2ccccc2Cl)CCCCC1=O"),
    ("N-nitroso-l-nebivolol", "OC(CN(N=O)CC(O)C1CCc2cc(F)ccc2O1)C1CCc2cc(F)ccc2O1"),
    ("N-nitrosomorpholine", "O=NN1CCOCC1"),
]


#  Single compound


class SingleTab(QWidget):
    def __init__(self, client: PubChemClient, limits: FDALimits, osr: OsrService):
        super().__init__()
        self._client = client
        self._limits = limits
        self._osr = osr
        self._worker: ScanWorker | None = None
        self._recogniser: RecogniseWorker | None = None
        self._result: ScanResult | None = None
        self.setAcceptDrops(True)

        # input row
        self._input = QLineEdit()
        self._input.setPlaceholderText(
            "SMILES, compound name, or CAS number.  Example: CCN(N=O)CC")
        self._input.setMinimumHeight(32)
        self._input.returnPressed.connect(self.assess)

        self._kind = QComboBox()
        self._kind.addItems(KINDS.keys())
        self._kind.setMinimumHeight(32)

        self._go = QPushButton("Assess")
        self._go.setMinimumHeight(32)
        self._go.setDefault(True)
        self._go.clicked.connect(self.assess)

        self._open_image = QPushButton("Open image…")
        self._open_image.setMinimumHeight(32)
        self._open_image.setToolTip(
            "Read the structure from a picture.\nYou can also drag an image onto "
            "this window.")
        self._open_image.clicked.connect(self._choose_image)
        if not osr_available():
            self._open_image.setEnabled(False)
            self._open_image.setToolTip(
                "Structure recognition is not installed.\n"
                "py -3.13 -m venv .venv-osr\n"
                ".venv-osr\\Scripts\\python -m pip install decimer")

        top = QHBoxLayout()
        top.addWidget(self._input, 1)
        top.addWidget(self._kind)
        top.addWidget(self._go)
        top.addWidget(self._open_image)

        examples = QHBoxLayout()
        examples.addWidget(QLabel("Try:"))
        for label, smiles in EXAMPLES:
            btn = QPushButton(label)
            btn.setFlat(True)
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
            btn.setStyleSheet("QPushButton { color: #1D4ED8; text-align: left; }")
            btn.clicked.connect(lambda _, s=smiles: self._use_example(s))
            examples.addWidget(btn)
        examples.addStretch(1)

        self._rules = QComboBox()
        self._rules.addItems(RULESETS.keys())
        self._rules.setToolTip(RULESET_HELP)
        self._rules.currentIndexChanged.connect(self._rules_changed)
        rules_label = QLabel("Rules:")
        rules_label.setToolTip(RULESET_HELP)
        examples.addWidget(rules_label)
        examples.addWidget(self._rules)

        # left: structure + dose
        self._structure = StructureView()

        self._mdd = QDoubleSpinBox()
        self._mdd.setRange(0, 100_000)
        self._mdd.setDecimals(1)
        self._mdd.setSingleStep(50)
        self._mdd.setValue(0)
        self._mdd.setSuffix("  mg/day")
        self._mdd.setSpecialValueText("not set")
        self._mdd.setToolTip(
            "Maximum daily dose of the drug product.\n"
            "Control limit (ppm) = AI (ng/day) / MDD (mg).")
        self._mdd.valueChanged.connect(self._refresh_verdict)

        dose = QHBoxLayout()
        dose.addWidget(QLabel("Max daily dose:"))
        dose.addWidget(self._mdd, 1)

        # Shown only after a picture is read, so a misreading can be corrected by hand
        # rather than silently producing a wrong number.
        self._read_smiles = QLineEdit()
        self._read_smiles.setPlaceholderText("structure read from the image")
        self._read_smiles.setToolTip(
            "What was read from the picture. If it is wrong, correct it here and "
            "press Re-score.")
        self._rescore = QPushButton("Re-score")
        self._rescore.clicked.connect(self._rescore_edited)
        self._read_row = QWidget()
        read_layout = QHBoxLayout(self._read_row)
        read_layout.setContentsMargins(0, 0, 0, 0)
        read_layout.addWidget(QLabel("Read as:"))
        read_layout.addWidget(self._read_smiles, 1)
        read_layout.addWidget(self._rescore)
        self._read_row.setVisible(False)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(self._structure)
        left_layout.addWidget(self._read_row)
        left_layout.addLayout(dose)
        left_layout.addStretch(1)

        # right: verdict + factors
        self._verdict = VerdictCard()
        self._factors = FactorTable()

        self._notes = QTextEdit()
        self._notes.setReadOnly(True)
        self._notes.setMaximumHeight(92)
        self._notes.setVisible(False)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.addWidget(self._verdict)
        right_layout.addWidget(QLabel("Factors in the calculation"))
        right_layout.addWidget(self._factors, 1)
        right_layout.addWidget(self._notes)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.addWidget(left)
        split.addWidget(right)
        split.setSizes([460, 640])

        # actions
        self._copy = QPushButton("Copy report")
        self._save = QPushButton("Save report…")
        for b in (self._copy, self._save):
            b.setEnabled(False)
        self._copy.clicked.connect(self._copy_report)
        self._save.clicked.connect(self._save_report)

        actions = QHBoxLayout()
        actions.addStretch(1)
        actions.addWidget(self._copy)
        actions.addWidget(self._save)

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addLayout(examples)
        layout.addWidget(split, 1)
        layout.addLayout(actions)

    # behaviour

    def _use_example(self, smiles: str) -> None:
        self._input.setText(smiles)
        self.assess()

    # reading a structure from a picture

    def dragEnterEvent(self, event) -> None:  # noqa: N802
        urls = event.mimeData().urls() if event.mimeData().hasUrls() else []
        if len(urls) == 1 and Path(urls[0].toLocalFile()).suffix.lower() in IMAGE_SUFFIXES:
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:  # noqa: N802
        path = Path(event.mimeData().urls()[0].toLocalFile())
        event.acceptProposedAction()
        self._recognise(path)

    def _choose_image(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open a structure image", str(Path.home()),
            "Images (*.png *.jpg *.jpeg *.bmp *.gif *.tif *.tiff);;All files (*)")
        if path:
            self._recognise(Path(path))

    def _recognise(self, path: Path) -> None:
        if self._recogniser and self._recogniser.isRunning():
            return
        if not osr_available():
            QMessageBox.information(self, "Not installed", self._open_image.toolTip())
            return

        self._input.setText(path.name)
        self._verdict.clear()
        self._verdict.show_error("Reading the structure from the picture...")
        self._factors.setRowCount(0)
        self._read_row.setVisible(False)
        self._structure.show_source_image(path)
        self._open_image.setEnabled(False)
        self._go.setEnabled(False)

        self._recogniser = RecogniseWorker(path, self._osr)
        self._recogniser.status.connect(self._verdict.show_error)
        self._recogniser.done.connect(self._on_recognised)
        self._recogniser.failed.connect(self._fail)
        self._recogniser.finished.connect(self._image_ready)
        self._recogniser.start()

    def _image_ready(self) -> None:
        self._open_image.setEnabled(osr_available())
        self._go.setEnabled(True)

    def _on_recognised(self, recognition) -> None:
        if not recognition.ok:
            self._verdict.show_error(
                recognition.error or "Nothing could be read from that picture.")
            return
        self._read_smiles.setText(recognition.smiles)
        self._read_row.setVisible(True)
        self._score_smiles(recognition.smiles)

    def _rescore_edited(self) -> None:
        text = self._read_smiles.text().strip()
        if text:
            self._score_smiles(text)

    def _ruleset(self) -> str:
        return RULESETS[self._rules.currentText()]

    def _rules_changed(self) -> None:
        """Re-score whatever is on screen, so switching rules shows its effect at once
        instead of waiting for the next lookup. Keeps the source picture on screen if
        the structure came from one."""
        if self._worker and self._worker.isRunning():
            return
        if not (self._result and self._result.smiles):
            return
        from_image = self._read_row.isVisible()
        self._score_smiles(self._result.smiles, keep_source=from_image)

    def _score_smiles(self, smiles: str, keep_source: bool = True) -> None:
        """Score a structure we already have, without going back to the network."""
        self._worker = ScanWorker(smiles, "smiles", self._client, self._limits,
                                  self._ruleset())
        self._worker.done.connect(
            self._show_keeping_source if keep_source else self._show)
        self._worker.failed.connect(self._fail)
        self._worker.start()

    def _show_keeping_source(self, result: ScanResult) -> None:
        """Same as _show, but leaves the source picture on screen for comparison."""
        self._result = result
        self._structure.show_smiles(result.smiles)
        self._factors.show_result(result)
        self._refresh_verdict()
        notes = list(result.notes) + (result.cpca.warnings if result.cpca else [])
        notes.append("Structure read from a picture by optical recognition. "
                     "Compare the two drawings before relying on this.")
        self._notes.setPlainText("\n".join(f"• {n}" for n in notes))
        self._notes.setVisible(True)
        self._copy.setEnabled(result.error is None)
        self._save.setEnabled(result.error is None)

    def assess(self) -> None:
        query = self._input.text().strip()
        if not query or (self._worker and self._worker.isRunning()):
            return

        self._go.setEnabled(False)
        self._go.setText("Working…")
        self._verdict.clear()
        self._factors.setRowCount(0)
        self._structure.clear()
        self._read_row.setVisible(False)
        self._notes.setVisible(False)

        self._worker = ScanWorker(query, KINDS[self._kind.currentText()],
                                  self._client, self._limits, self._ruleset())
        self._worker.done.connect(self._show)
        self._worker.failed.connect(self._fail)
        self._worker.finished.connect(self._ready)
        self._worker.start()

    def _ready(self) -> None:
        self._go.setEnabled(True)
        self._go.setText("Assess")

    def _fail(self, message: str) -> None:
        self._verdict.show_error(message)
        self._copy.setEnabled(False)
        self._save.setEnabled(False)

    def _show(self, result: ScanResult) -> None:
        self._result = result
        self._structure.show_smiles(result.smiles)
        self._factors.show_result(result)
        self._refresh_verdict()

        notes = list(result.notes)
        if result.cpca:
            notes += result.cpca.warnings
        if result.resolved_name:
            notes.append(f"Resolved to: {result.resolved_name}")
        if result.smiles and result.input_kind != "smiles":
            notes.append(f"SMILES: {result.smiles}")
        self._notes.setPlainText("\n".join(f"• {n}" for n in notes))
        self._notes.setVisible(bool(notes))

        ok = result.error is None
        self._copy.setEnabled(ok)
        self._save.setEnabled(ok)

    def _refresh_verdict(self) -> None:
        if self._result is not None:
            self._verdict.show(self._result, self._mdd.value() or None)

    def _report_text(self) -> str:
        return format_scan_report(self._result, self._mdd.value() or None)

    def _copy_report(self) -> None:
        QGuiApplication.clipboard().setText(self._report_text())
        self._go.setText("Copied")
        self._go.setText("Assess")

    def _save_report(self) -> None:
        name = (self._result.query[:40] if self._result else "report")
        safe = "".join(c if c.isalnum() or c in "-_ " else "_" for c in name).strip()
        path, _ = QFileDialog.getSaveFileName(
            self, "Save report", str(Path.home() / f"{safe or 'report'}.txt"),
            "Text files (*.txt);;All files (*)")
        if path:
            Path(path).write_text(self._report_text(), encoding="utf-8")


#  Batch


class BatchTab(QWidget):
    COLUMNS = ["Input", "Category", "AI (ng/day)", "Basis", "Score", "Notes"]

    def __init__(self) -> None:
        super().__init__()
        self._worker: BatchWorker | None = None
        self._results: list[ScanResult | None] = []
        self._queries: list[str] = []

        self._open = QPushButton("Open file…")
        self._open.clicked.connect(self._open_file)

        self._kind = QComboBox()
        self._kind.addItems(KINDS.keys())

        self._run = QPushButton("Assess all")
        self._run.setEnabled(False)
        self._run.clicked.connect(self._run_batch)

        self._cancel = QPushButton("Cancel")
        self._cancel.setEnabled(False)
        self._cancel.clicked.connect(self._cancel_batch)

        self._export = QPushButton("Export CSV…")
        self._export.setEnabled(False)
        self._export.clicked.connect(self._export_csv)

        self._rules = QComboBox()
        self._rules.addItems(RULESETS.keys())
        self._rules.setToolTip(RULESET_HELP)
        rules_label = QLabel("Rules:")
        rules_label.setToolTip(RULESET_HELP)

        top = QHBoxLayout()
        top.addWidget(self._open)
        top.addWidget(QLabel("Treat entries as:"))
        top.addWidget(self._kind)
        top.addWidget(rules_label)
        top.addWidget(self._rules)
        top.addWidget(self._run)
        top.addWidget(self._cancel)
        top.addStretch(1)
        top.addWidget(self._export)

        self._status = QLabel(
            "Open a .csv, .xlsx or .txt file. The first column with usable values is "
            "read; one compound per row.")
        self._status.setWordWrap(True)
        self._status.setStyleSheet("color: #57534E;")

        self._progress = QProgressBar()
        self._progress.setVisible(False)

        self._table = QTableWidget(0, len(self.COLUMNS))
        self._table.setHorizontalHeaderLabels(self.COLUMNS)
        self._table.verticalHeader().setVisible(False)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setAlternatingRowColors(True)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.Stretch)

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self._status)
        layout.addWidget(self._progress)
        layout.addWidget(self._table, 1)

    # loading

    def _open_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open compound list", str(Path.home()),
            "Compound lists (*.csv *.txt *.xlsx *.xlsm);;All files (*)")
        if not path:
            return
        try:
            self._queries = read_compound_list(Path(path))
        except Exception as e:  # noqa: BLE001
            QMessageBox.warning(self, "Could not read file", f"{type(e).__name__}: {e}")
            return

        if not self._queries:
            QMessageBox.information(self, "Nothing found",
                                    "No usable values were found in that file.")
            return

        self._results = [None] * len(self._queries)
        self._table.setRowCount(len(self._queries))
        for i, q in enumerate(self._queries):
            self._table.setItem(i, 0, QTableWidgetItem(q))
            for col in range(1, len(self.COLUMNS)):
                self._table.setItem(i, col, QTableWidgetItem(""))
        self._status.setText(
            f"{len(self._queries)} compounds loaded from {Path(path).name}")
        self._run.setEnabled(True)
        self._export.setEnabled(False)

    # running

    def _run_batch(self) -> None:
        self._run.setEnabled(False)
        self._cancel.setEnabled(True)
        self._progress.setVisible(True)
        self._progress.setRange(0, len(self._queries))
        self._progress.setValue(0)

        self._worker = BatchWorker(self._queries, KINDS[self._kind.currentText()],
                                   RULESETS[self._rules.currentText()])
        self._worker.progress.connect(self._on_progress)
        self._worker.row_done.connect(self._on_row)
        self._worker.finished_all.connect(self._on_finished)
        self._worker.failed.connect(
            lambda m: QMessageBox.warning(self, "Batch failed", m))
        self._worker.start()

    def _cancel_batch(self) -> None:
        if self._worker:
            self._worker.cancel()
        self._cancel.setEnabled(False)

    def _on_progress(self, done: int, total: int, label: str) -> None:
        self._progress.setValue(done)
        self._status.setText(f"{done}/{total}  ·  {label}")

    def _on_row(self, index: int, result: ScanResult) -> None:
        self._results[index] = result
        site = result.cpca.governing_site if result.cpca else None
        values = [
            str(result.final_category or ""),
            f"{result.final_ai_ng_day:g}" if result.final_ai_ng_day else "",
            result.final_basis if not result.error else "",
            str(site.potency_score) if site and site.potency_score is not None else "",
            result.error or "; ".join(result.notes)[:120],
        ]
        for offset, value in enumerate(values, start=1):
            self._table.setItem(index, offset, QTableWidgetItem(value))
        self._progress.setValue(index + 1)

    def _on_finished(self, completed: int) -> None:
        self._progress.setVisible(False)
        self._cancel.setEnabled(False)
        self._run.setEnabled(True)
        self._export.setEnabled(True)
        scored = sum(1 for r in self._results if r and r.final_ai_ng_day)
        strict = sum(1 for r in self._results
                     if r and r.final_ai_ng_day and r.final_ai_ng_day <= 100)
        self._status.setText(
            f"Done. {completed} assessed, {scored} with a limit, "
            f"{strict} at 100 ng/day or tighter.")

    def _export_csv(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Export results", str(Path.home() / "impurity_results.csv"),
            "CSV files (*.csv)")
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            writer = csv.writer(fh)
            writer.writerow(["input", "smiles", "category", "ai_ng_day", "basis",
                             "potency_score", "factors", "notes"])
            for query, r in zip(self._queries, self._results):
                if r is None:
                    writer.writerow([query, "", "", "", "not assessed", "", "", ""])
                    continue
                site = r.cpca.governing_site if r.cpca else None
                factors = "; ".join(f"{f.score:+d} {f.name}"
                                    for f in (site.factors if site else []))
                writer.writerow([
                    query, r.smiles or "", r.final_category or "",
                    r.final_ai_ng_day or "", r.final_basis,
                    site.potency_score if site else "",
                    factors or (site.triage_exit if site else ""),
                    r.error or "; ".join(r.notes),
                ])


def read_compound_list(path: Path) -> list[str]:
    """Read one compound per row from csv / txt / xlsx, using the first usable column."""
    suffix = path.suffix.lower()

    if suffix in (".xlsx", ".xlsm"):
        import openpyxl
        ws = openpyxl.load_workbook(path, data_only=True).worksheets[0]
        rows = [[("" if c.value is None else str(c.value).strip()) for c in row]
                for row in ws.iter_rows()]
    elif suffix == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as fh:
            rows = [[c.strip() for c in row] for row in csv.reader(fh)]
    else:
        rows = [[line.strip()] for line in
                path.read_text(encoding="utf-8", errors="replace").splitlines()]

    if not rows:
        return []

    width = max(len(r) for r in rows)
    best_col, best_count = 0, -1
    for col in range(width):
        values = [r[col] for r in rows if col < len(r) and r[col]]
        # Prefer the column with the most non-empty, non-numeric entries.
        usable = [v for v in values if not v.replace(".", "", 1).isdigit()]
        if len(usable) > best_count:
            best_col, best_count = col, len(usable)

    out, seen = [], set()
    for i, row in enumerate(rows):
        if best_col >= len(row):
            continue
        value = row[best_col]
        if not value or value.lower() in seen:
            continue
        if i == 0 and value.lower() in {"smiles", "name", "compound", "cas",
                                        "structure", "impurity", "impurity_name"}:
            continue  # header
        seen.add(value.lower())
        out.append(value)
    return out


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Nitrosamine Impurity Scanner (FDA CPCA)")
        self.resize(1180, 780)

        client = PubChemClient(verbose=False)
        limits = FDALimits()
        self._osr = OsrService()

        tabs = QTabWidget()
        tabs.addTab(SingleTab(client, limits, self._osr), "Single compound")
        tabs.addTab(BatchTab(), "Batch")
        self.setCentralWidget(tabs)

        bar = self.statusBar()
        stale = limits.staleness_warning()
        recognition = ("image recognition ready" if osr_available()
                       else "image recognition not installed")
        bar.showMessage(
            (stale or f"FDA published limits: {len(limits)} entries, "
                      f"snapshot {limits.fetched}")
            + f"   ·   {recognition}"
            + "   ·   Screening aid for a qualified assessor, not a regulatory determination.")

    def closeEvent(self, event) -> None:  # noqa: N802
        """Shut the recognition service down; it holds a loaded model in memory."""
        self._osr.stop()
        super().closeEvent(event)


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("Nitrosamine Impurity Scanner")
    app.setFont(QFont("Segoe UI", 10))
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
