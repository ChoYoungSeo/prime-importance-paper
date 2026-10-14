#!/usr/bin/env python3
"""Build the Concrete p=8 figure from the pinned 100-subsample summary.

Run from any working directory. Requires ReportLab and a supported TrueType
font pair. The release figure used Arial; other fonts preserve the plotted
values but change the PDF bytes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

from reportlab.lib.colors import Color, HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CSV = ROOT / "results/submission_20260920/concrete/distinct_824/coordinate_summary.csv"
DEFAULT_OUTPUT = ROOT / "manuscript/figures/concrete_p8_distinct_subsamples.pdf"
DEFAULT_RECEIPT = ROOT / "results/submission_20260920/concrete/FIGURE_RECEIPT.json"
EXPECTED_CSV_SHA256 = "be2d8c332feae206cf7242f8890b9290e7be6a2f8c199659b9c5b34f09bb0162"

METHODS = ["PRIME", "Full LL-PME", "PME-kNN", "FIRST-SHAP", "Shapley-kNN"]
VARIABLES = [
    "cement", "slag", "ash", "water_binder_ratio", "superplastic",
    "coarseagg", "fineagg", "age",
]
LABELS = [
    "Cement", "Slag", "Ash", "Water / binder ratio", "Superplasticizer",
    "Coarse aggregate", "Fine aggregate", "Age",
]
COLORS = {
    "PRIME": HexColor("#205A7B"),
    "Full LL-PME": HexColor("#557B75"),
    "PME-kNN": HexColor("#785F88"),
    "FIRST-SHAP": HexColor("#B95F3D"),
    "Shapley-kNN": HexColor("#785F88"),
}
SYMBOLS = {
    "PRIME": "circle", "Full LL-PME": "square", "PME-kNN": "diamond",
    "FIRST-SHAP": "circle", "Shapley-kNN": "square",
}
OFFSETS = {
    "PRIME": 3.65, "Full LL-PME": 0.0, "PME-kNN": -3.65,
    "FIRST-SHAP": 2.4, "Shapley-kNN": -2.4,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_data(path: Path, expected_n: int) -> dict[tuple[str, str], tuple[float, float, float]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    expected = {(method, variable) for method in METHODS for variable in VARIABLES}
    if len(rows) != 40:
        raise ValueError(f"Expected 40 method-coordinate rows; got {len(rows)}")
    data = {}
    for row in rows:
        key = (row["method"], row["variable_name"])
        if key not in expected or key in data:
            raise ValueError(f"Unexpected or duplicate method-coordinate row: {key}")
        observed_n = row.get("n") or row.get("n_complete")
        if observed_n is None or int(observed_n) != expected_n:
            raise ValueError(f"Expected {expected_n} complete resamples for {key}")
        values = tuple(float(row[name]) for name in ("q05", "median", "q95"))
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"Nonfinite plotted value for {key}")
        if not 0 <= values[0] <= values[1] <= values[2] <= 0.8:
            raise ValueError(f"Invalid 5th/median/95th percentile order for {key}")
        data[key] = values
    if set(data) != expected:
        raise ValueError("Missing method-coordinate rows")
    return data


def select_fonts() -> tuple[Path, Path]:
    candidates = [
        (
            Path("/System/Library/Fonts/Supplemental/Arial.ttf"),
            Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf"),
        ),
        (
            Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"),
            Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
        ),
        (
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        ),
    ]
    for regular, bold in candidates:
        if regular.is_file() and bold.is_file():
            return regular, bold
    raise FileNotFoundError("No supported regular/bold TrueType font pair found")


def draw_figure(
    output: Path,
    data: dict[tuple[str, str], tuple[float, float, float]],
    source_sha256: str,
) -> tuple[Path, Path]:
    regular, bold = select_fonts()
    pdfmetrics.registerFont(TTFont("FigureRegular", str(regular)))
    pdfmetrics.registerFont(TTFont("FigureBold", str(bold)))

    output.parent.mkdir(parents=True, exist_ok=True)
    width, height = 468.0, 270.0  # 6.5 x 3.75 in; width equals the manuscript text width.
    page = canvas.Canvas(str(output), pagesize=(width, height), pageCompression=1, invariant=1)
    page.setTitle("Concrete p=8: five-method coordinate allocations")
    page.setSubject(
        "Median and empirical 5th–95th percentiles; source coordinate_summary.csv "
        "SHA-256 " + source_sha256
    )
    page.setAuthor("NJI")
    page.setCreator("ReportLab")

    ink = HexColor("#28323B")
    muted = HexColor("#53606A")
    axis = HexColor("#9AA8B3")
    x0, x1, x2, x3 = 102.0, 268.0, 292.0, 458.0
    centers = [203.0 - 19.0 * index for index in range(8)]

    def xcoord(panel: int, value: float) -> float:
        left, right = (x0, x1) if panel == 0 else (x2, x3)
        return left + value / 0.8 * (right - left)

    def line(xa: float, ya: float, xb: float, yb: float, color: Color, stroke: float = 0.5) -> None:
        page.setStrokeColor(color)
        page.setLineWidth(stroke)
        page.line(xa, ya, xb, yb)

    def marker(x: float, y: float, color: Color, shape: str, size: float = 2.75) -> None:
        page.setStrokeColor(color)
        page.setFillColor(Color(1, 1, 1))
        page.setLineWidth(1.2)
        if shape == "circle":
            page.circle(x, y, size, stroke=1, fill=1)
        elif shape == "square":
            page.rect(x - size, y - size, 2 * size, 2 * size, stroke=1, fill=1)
        else:
            diamond = page.beginPath()
            diamond.moveTo(x, y + size + 0.2)
            diamond.lineTo(x + size + 0.2, y)
            diamond.lineTo(x, y - size - 0.2)
            diamond.lineTo(x - size - 0.2, y)
            diamond.close()
            page.drawPath(diamond, stroke=1, fill=1)

    def legend(x: float, y: float, method: str) -> None:
        color = COLORS[method]
        line(x, y, x + 12, y, color, 1.3)
        marker(x + 6, y, color, SYMBOLS[method], 2.45)
        page.setFont("FigureRegular", 7.25)
        page.setFillColor(ink)
        page.drawString(x + 17, y - 2.4, method)

    page.setFillColor(ink)
    page.setFont("FigureBold", 9.0)
    page.drawString(x0, 253.0, "PME-target methods")
    page.drawString(x2, 253.0, "Shapley-target methods")
    legend(x0, 236.0, "PRIME")
    legend(x0 + 65, 236.0, "Full LL-PME")
    legend(x0, 224.0, "PME-kNN")
    legend(x2, 236.0, "FIRST-SHAP")
    legend(x2 + 79, 236.0, "Shapley-kNN")

    for index, y in enumerate(centers):
        page.setFillColor(ink)
        page.setFont("FigureRegular", 7.8)
        page.drawRightString(94.0, y - 2.6, LABELS[index])
    for panel, left, right in ((0, x0, x1), (1, x2, x3)):
        for tick in (0, 0.2, 0.4, 0.6, 0.8):
            x = xcoord(panel, tick)
            if tick == 0:
                line(x, 61.0, x, 212.0, axis, 0.65)
            else:
                line(x, 61.0, x, 64.0, axis, 0.65)
            page.setFillColor(muted)
            page.setFont("FigureRegular", 7.2)
            page.drawCentredString(x, 48.5, f"{tick:.1f}")
        line(left, 61.0, right, 61.0, axis, 0.65)
        page.setFillColor(muted)
        page.setFont("FigureRegular", 7.4)
        page.drawCentredString((left + right) / 2, 31.5, "Relative allocation")

    for panel, group in ((0, METHODS[:3]), (1, METHODS[3:])):
        for index, variable in enumerate(VARIABLES):
            for method in group:
                q05, median, q95 = data[(method, variable)]
                y = centers[index] + OFFSETS[method]
                xa, xm, xb = (xcoord(panel, value) for value in (q05, median, q95))
                color = COLORS[method]
                line(xa, y, xb, y, color, 1.15)
                line(xa, y - 2.25, xa, y + 2.25, color, 0.85)
                line(xb, y - 2.25, xb, y + 2.25, color, 0.85)
                marker(xm, y, color, SYMBOLS[method], 2.65)

    page.showPage()
    page.save()
    return regular, bold


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV, help="Validated coordinate_summary.csv")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output manuscript PDF")
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT, help="JSON hash receipt")
    parser.add_argument("--n", type=int, default=100, help="Complete paired resamples per coordinate")
    parser.add_argument(
        "--design", default="100 paired 824/1030 (80%) without-replacement subsamples",
        help="Design description recorded in the hash receipt",
    )
    parser.add_argument(
        "--source-sha256", default=EXPECTED_CSV_SHA256,
        help="Required SHA-256 of the validated input CSV",
    )
    args = parser.parse_args()

    if args.n < 1:
        parser.error("--n must be positive")
    csv_path = args.csv.resolve()
    output = args.output.resolve()
    receipt_path = args.receipt.resolve()
    source_hash = sha256(csv_path)
    if source_hash != args.source_sha256.lower():
        raise ValueError(
            "CSV hash differs from the specified validated source: " + source_hash
        )
    data = load_data(csv_path, args.n)
    regular, bold = draw_figure(output, data, source_hash)
    receipt = {
        "builder": Path(__file__).name,
        "source_csv": os.path.relpath(csv_path, ROOT),
        "source_sha256": source_hash,
        "replications": args.n,
        "design": args.design,
        "output_pdf": os.path.relpath(output, ROOT),
        "output_sha256": sha256(output),
        "regular_font": regular.name,
        "bold_font": bold.name,
        "size_pt": [468, 270],
        "median_markers": len(data),
        "percentile_ranges": len(data),
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Receipt {receipt_path}")
    print(f"Source SHA-256 {source_hash}")
    print(f"Output SHA-256 {receipt['output_sha256']}")


if __name__ == "__main__":
    main()
