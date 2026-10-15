#!/usr/bin/env python3
"""Render a targeted rooted tree as a compact vector PDF."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


class SimplePDFCanvas:
    """Small dependency-free PDF writer for lines and built-in Helvetica text."""

    def __init__(self, output: Path, width: float, height: float) -> None:
        self.output = output
        self.width = width
        self.height = height
        self.commands: list[str] = []
        self.font = "F1"
        self.font_size = 8.0
        self.title = ""

    @staticmethod
    def _number(value: float) -> str:
        return f"{value:.3f}".rstrip("0").rstrip(".")

    @staticmethod
    def _text(value: str) -> str:
        safe = value.encode("latin-1", errors="replace").decode("latin-1")
        return safe.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")

    def set_title(self, title: str) -> None:
        self.title = title

    def set_font(self, name: str, size: float) -> None:
        self.font = "F2" if "Bold" in name else "F1"
        self.font_size = size

    def set_stroke_color(self, red: float, green: float, blue: float) -> None:
        self.commands.append(
            f"{self._number(red)} {self._number(green)} {self._number(blue)} RG"
        )

    def set_fill_color(self, red: float, green: float, blue: float) -> None:
        self.commands.append(
            f"{self._number(red)} {self._number(green)} {self._number(blue)} rg"
        )

    def set_line_width(self, width: float) -> None:
        self.commands.append(f"{self._number(width)} w")

    def line(self, x1: float, y1: float, x2: float, y2: float) -> None:
        self.commands.append(
            f"{self._number(x1)} {self._number(y1)} m "
            f"{self._number(x2)} {self._number(y2)} l S"
        )

    def draw_string(self, x: float, y: float, value: str) -> None:
        self.commands.append(
            f"BT /{self.font} {self._number(self.font_size)} Tf "
            f"{self._number(x)} {self._number(y)} Td ({self._text(value)}) Tj ET"
        )

    def _estimated_width(self, value: str) -> float:
        return len(value) * self.font_size * 0.52

    def draw_right_string(self, x: float, y: float, value: str) -> None:
        self.draw_string(x - self._estimated_width(value), y, value)

    def draw_centred_string(self, x: float, y: float, value: str) -> None:
        self.draw_string(x - self._estimated_width(value) / 2, y, value)

    def save(self) -> None:
        stream = ("\n".join(self.commands) + "\n").encode("latin-1", errors="replace")
        title = self._text(self.title)
        objects = [
            b"<< /Type /Catalog /Pages 2 0 R >>",
            b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
            (
                f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self._number(self.width)} "
                f"{self._number(self.height)}] /Resources << /Font << /F1 5 0 R /F2 6 0 R >> >> "
                "/Contents 4 0 R >>"
            ).encode("ascii"),
            b"<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n" + stream + b"endstream",
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>",
            f"<< /Title ({title}) >>".encode("latin-1", errors="replace"),
        ]
        document = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = [0]
        for number, obj in enumerate(objects, start=1):
            offsets.append(len(document))
            document.extend(f"{number} 0 obj\n".encode("ascii"))
            document.extend(obj)
            document.extend(b"\nendobj\n")
        xref = len(document)
        document.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
        document.extend(b"0000000000 65535 f \n")
        for offset in offsets[1:]:
            document.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
        document.extend(
            (
                f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R /Info 7 0 R >>\n"
                f"startxref\n{xref}\n%%EOF\n"
            ).encode("ascii")
        )
        self.output.write_bytes(document)


COLORS = {
    "BGN": (0.12, 0.40, 0.75),
    "ASPN": (0.80, 0.25, 0.20),
    "DCN": (0.05, 0.60, 0.60),
    "DCN2": (0.08, 0.38, 0.42),
    "FMOD": (0.40, 0.55, 0.18),
    "FMODA": (0.55, 0.68, 0.25),
    "FMODB": (0.27, 0.48, 0.12),
    "FMODB_CHECK": (0.15, 0.30, 0.05),
    "LUM": (0.65, 0.35, 0.75),
    "LUM_CHECK": (0.42, 0.15, 0.55),
    "LUML": (0.82, 0.55, 0.88),
    "PRELP": (0.75, 0.50, 0.05),
    "KERA": (0.55, 0.35, 0.10),
    "OMD": (0.75, 0.65, 0.20),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", choices=("class1", "class2"), required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    nodes = read_tsv(args.nodes)
    edges = read_tsv(args.edges)
    tips = [node for node in nodes if node["node_type"] == "tip"]
    maximum_x = max(float(node["x"]) for node in nodes)
    maximum_y = max(float(node["y"]) for node in nodes)
    row_height = 12.5
    top = 85.0
    bottom = 45.0
    width = 1500.0
    height = max(800.0, top + bottom + row_height * max(1, len(tips) - 1))
    tree_left = 45.0
    tree_right = 880.0
    label_left = tree_right + 12.0

    def px(value: str) -> float:
        return tree_left + (float(value) / maximum_x) * (tree_right - tree_left)

    def py(value: str) -> float:
        return height - top - (float(value) / maximum_y) * row_height * maximum_y

    title = (
        "Expanded BGN-DCN family tree"
        if args.analysis == "class1"
        else "Expanded FMOD-LUM family tree"
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pdf = SimplePDFCanvas(args.output, width, height)
    pdf.set_title(title)
    pdf.set_font("Helvetica-Bold", 15)
    pdf.draw_string(tree_left, height - 30, title)
    pdf.set_font("Helvetica", 8)
    pdf.draw_string(tree_left, height - 47, "Branch labels show SH-aLRT/UFBoot support when both values are available.")

    pdf.set_stroke_color(0.55, 0.55, 0.55)
    pdf.set_line_width(0.45)
    for edge in edges:
        x1, y1 = px(edge["parent_x"]), py(edge["parent_y"])
        x2, y2 = px(edge["child_x"]), py(edge["child_y"])
        pdf.line(x1, y1, x1, y2)
        pdf.line(x1, y2, x2, y2)

    for node in nodes:
        x, y = px(node["x"]), py(node["y"])
        if node["node_type"] == "tip":
            group = node["group"]
            color = COLORS.get(group, (0.25, 0.25, 0.25))
            pdf.set_fill_color(*color)
            parts = node["label"].split("|")
            if len(parts) >= 3:
                parts[1] = parts[1].replace("_", " ")
            label = "|".join(parts)
            pdf.set_font("Helvetica-Bold" if group.endswith("_CHECK") else "Helvetica", 10.0)
            pdf.draw_string(label_left, y - 2.8, label)
        else:
            support = node["support"]
            if "/" not in support:
                continue
            try:
                alrt, ufboot = (float(value) for value in support.split("/", 1))
            except ValueError:
                continue
            if alrt >= 80 and ufboot >= 90:
                pdf.set_fill_color(0.15, 0.15, 0.15)
                pdf.set_font("Helvetica", 5.4)
                pdf.draw_right_string(x - 2, y + 1.5, support)

    pdf.set_stroke_color(0, 0, 0)
    scale_length = 0.1
    scale_pixels = (scale_length / maximum_x) * (tree_right - tree_left)
    pdf.line(tree_left, 24, tree_left + scale_pixels, 24)
    pdf.set_font("Helvetica", 7)
    pdf.draw_centred_string(tree_left + scale_pixels / 2, 12, "0.1 substitutions/site")
    pdf.save()


if __name__ == "__main__":
    main()
