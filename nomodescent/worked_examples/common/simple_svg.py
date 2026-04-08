from __future__ import annotations

from pathlib import Path

import numpy as np


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def line_chart_svg(
    path: Path,
    x: np.ndarray,
    series: list[tuple[str, np.ndarray, str]],
    *,
    width: int = 860,
    height: int = 520,
    title: str = "",
    x_label: str = "",
    y_label: str = "",
) -> None:
    margin_left = 80
    margin_right = 24
    margin_top = 50
    margin_bottom = 70
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    y_values = np.concatenate([values for _, values, _ in series])
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    y_min = float(np.min(y_values))
    y_max = float(np.max(y_values))
    if abs(x_max - x_min) < 1e-15:
        x_max = x_min + 1.0
    if abs(y_max - y_min) < 1e-15:
        y_max = y_min + 1.0
    y_pad = 0.05 * (y_max - y_min)
    y_min -= y_pad
    y_max += y_pad

    def xp(value: float) -> float:
        return margin_left + (value - x_min) / (x_max - x_min) * plot_width

    def yp(value: float) -> float:
        return margin_top + plot_height - (value - y_min) / (y_max - y_min) * plot_height

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width/2:.1f}" y="28" text-anchor="middle" font-size="20" font-family="monospace">{title}</text>',
        f'<line x1="{margin_left}" y1="{margin_top + plot_height}" x2="{margin_left + plot_width}" y2="{margin_top + plot_height}" stroke="#222" stroke-width="1.5"/>',
        f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_height}" stroke="#222" stroke-width="1.5"/>',
    ]
    for label, values, color in series:
        points = " ".join(f"{xp(float(xv)):.2f},{yp(float(yv)):.2f}" for xv, yv in zip(x, values, strict=True))
        elements.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.5" points="{points}"/>')
        for xv, yv in zip(x, values, strict=True):
            elements.append(f'<circle cx="{xp(float(xv)):.2f}" cy="{yp(float(yv)):.2f}" r="2.2" fill="{color}"/>')
    elements.append(f'<text x="{margin_left + plot_width/2:.1f}" y="{height - 16}" text-anchor="middle" font-size="14" font-family="monospace">{x_label}</text>')
    elements.append(f'<text x="18" y="{margin_top + plot_height/2:.1f}" text-anchor="middle" transform="rotate(-90 18 {margin_top + plot_height/2:.1f})" font-size="14" font-family="monospace">{y_label}</text>')
    elements.append("</svg>")
    _write(path, "\n".join(elements))


def bar_chart_svg(
    path: Path,
    labels: list[str],
    values: np.ndarray,
    *,
    width: int = 860,
    height: int = 520,
    title: str = "",
    y_label: str = "",
    color: str = "#4477aa",
) -> None:
    margin_left = 80
    margin_right = 24
    margin_top = 50
    margin_bottom = 80
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    ymax = max(1.0, float(np.max(values)) * 1.1)
    bar_w = plot_width / max(1, len(values))
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width/2:.1f}" y="28" text-anchor="middle" font-size="20" font-family="monospace">{title}</text>',
    ]
    for idx, (label, value) in enumerate(zip(labels, values, strict=True)):
        x0 = margin_left + idx * bar_w + 0.18 * bar_w
        bw = 0.64 * bar_w
        bh = plot_height * float(value) / ymax
        y0 = margin_top + plot_height - bh
        elements.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{bw:.2f}" height="{bh:.2f}" fill="{color}" stroke="#355c7d"/>')
        elements.append(f'<text x="{x0 + bw/2:.2f}" y="{margin_top + plot_height + 24}" text-anchor="middle" font-size="12" font-family="monospace">{label}</text>')
        elements.append(f'<text x="{x0 + bw/2:.2f}" y="{y0 - 6:.2f}" text-anchor="middle" font-size="12" font-family="monospace">{float(value):.3g}</text>')
    elements.append(f'<text x="18" y="{margin_top + plot_height/2:.1f}" text-anchor="middle" transform="rotate(-90 18 {margin_top + plot_height/2:.1f})" font-size="14" font-family="monospace">{y_label}</text>')
    elements.append("</svg>")
    _write(path, "\n".join(elements))
