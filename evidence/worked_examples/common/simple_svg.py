from __future__ import annotations

from pathlib import Path

import numpy as np


def bar_chart_svg(
    path: Path,
    labels: list[str],
    values: np.ndarray,
    *,
    title: str,
    y_label: str,
    color: str = "#355070",
) -> None:
    values = np.asarray(values, dtype=float)
    width = 480
    height = 320
    margin_left = 70
    margin_bottom = 55
    margin_top = 40
    chart_width = width - margin_left - 20
    chart_height = height - margin_top - margin_bottom
    vmax = float(values.max()) if values.size else 1.0
    vmax = 1.0 if vmax <= 0.0 else vmax
    bar_width = chart_width / max(len(labels), 1)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        f'<text x="{width / 2}" y="24" text-anchor="middle" font-size="16">{title}</text>',
        f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height - margin_bottom}" stroke="#222"/>',
        f'<line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - 20}" y2="{height - margin_bottom}" stroke="#222"/>',
        f'<text x="20" y="{margin_top + chart_height / 2}" transform="rotate(-90 20,{margin_top + chart_height / 2})" text-anchor="middle" font-size="12">{y_label}</text>',
    ]

    for idx, (label, value) in enumerate(zip(labels, values, strict=True)):
        x = margin_left + idx * bar_width + 12
        bar_h = chart_height * (value / vmax)
        y = margin_top + chart_height - bar_h
        parts.append(f'<rect x="{x}" y="{y}" width="{bar_width - 24}" height="{bar_h}" fill="{color}"/>')
        parts.append(f'<text x="{x + (bar_width - 24) / 2}" y="{height - margin_bottom + 16}" text-anchor="middle" font-size="11">{label}</text>')
        parts.append(f'<text x="{x + (bar_width - 24) / 2}" y="{y - 6}" text-anchor="middle" font-size="11">{value:.3f}</text>')

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")
