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

    for i in range(5):
        xt = x_min + i * (x_max - x_min) / 4.0
        px = xp(xt)
        elements.append(f'<line x1="{px:.2f}" y1="{margin_top + plot_height}" x2="{px:.2f}" y2="{margin_top + plot_height + 6}" stroke="#222"/>')
        elements.append(f'<text x="{px:.2f}" y="{margin_top + plot_height + 24}" text-anchor="middle" font-size="12" font-family="monospace">{xt:.3g}</text>')
    for i in range(5):
        yt = y_min + i * (y_max - y_min) / 4.0
        py = yp(yt)
        elements.append(f'<line x1="{margin_left - 6}" y1="{py:.2f}" x2="{margin_left}" y2="{py:.2f}" stroke="#222"/>')
        elements.append(f'<text x="{margin_left - 10}" y="{py + 4:.2f}" text-anchor="end" font-size="12" font-family="monospace">{yt:.3g}</text>')
        elements.append(f'<line x1="{margin_left}" y1="{py:.2f}" x2="{margin_left + plot_width}" y2="{py:.2f}" stroke="#ddd" stroke-width="1"/>')

    for label, values, color in series:
        points = " ".join(f"{xp(float(xv)):.2f},{yp(float(yv)):.2f}" for xv, yv in zip(x, values, strict=True))
        elements.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.5" points="{points}"/>')
        for xv, yv in zip(x, values, strict=True):
            elements.append(f'<circle cx="{xp(float(xv)):.2f}" cy="{yp(float(yv)):.2f}" r="2.2" fill="{color}"/>')

    legend_x = width - 230
    legend_y = margin_top + 10
    elements.append(f'<rect x="{legend_x}" y="{legend_y}" width="190" height="{28 * len(series) + 10}" fill="#fafafa" stroke="#bbb"/>')
    for idx, (label, _values, color) in enumerate(series):
        y0 = legend_y + 22 + 26 * idx
        elements.append(f'<line x1="{legend_x + 14}" y1="{y0}" x2="{legend_x + 42}" y2="{y0}" stroke="{color}" stroke-width="2.5"/>')
        elements.append(f'<text x="{legend_x + 52}" y="{y0 + 4}" font-size="13" font-family="monospace">{label}</text>')

    elements.append(f'<text x="{margin_left + plot_width/2:.1f}" y="{height - 16}" text-anchor="middle" font-size="14" font-family="monospace">{x_label}</text>')
    elements.append(f'<text x="18" y="{margin_top + plot_height/2:.1f}" text-anchor="middle" transform="rotate(-90 18 {margin_top + plot_height/2:.1f})" font-size="14" font-family="monospace">{y_label}</text>')
    elements.append("</svg>")
    _write(path, "\n".join(elements))


def heatmap_svg(
    path: Path,
    x_values: np.ndarray,
    y_values: np.ndarray,
    grid: np.ndarray,
    palette: dict[int, str],
    labels: dict[int, str],
    *,
    width: int = 860,
    height: int = 560,
    title: str = "",
    x_label: str = "",
    y_label: str = "",
) -> None:
    margin_left = 80
    margin_right = 30
    margin_top = 50
    margin_bottom = 70
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    nx = len(x_values)
    ny = len(y_values)
    cell_w = plot_width / nx
    cell_h = plot_height / ny

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width/2:.1f}" y="28" text-anchor="middle" font-size="20" font-family="monospace">{title}</text>',
    ]

    for iy in range(ny):
        for ix in range(nx):
            cls = int(grid[iy, ix])
            color = palette.get(cls, "#ffffff")
            x0 = margin_left + ix * cell_w
            y0 = margin_top + (ny - 1 - iy) * cell_h
            elements.append(
                f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{cell_w + 0.2:.2f}" height="{cell_h + 0.2:.2f}" fill="{color}" stroke="#ffffff" stroke-width="0.4"/>'
            )

    elements.append(f'<rect x="{margin_left}" y="{margin_top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#222" stroke-width="1.4"/>')
    for i in range(5):
        xt = float(x_values[0] + i * (x_values[-1] - x_values[0]) / 4.0)
        px = margin_left + i * plot_width / 4.0
        elements.append(f'<line x1="{px:.2f}" y1="{margin_top + plot_height}" x2="{px:.2f}" y2="{margin_top + plot_height + 6}" stroke="#222"/>')
        elements.append(f'<text x="{px:.2f}" y="{margin_top + plot_height + 24}" text-anchor="middle" font-size="12" font-family="monospace">{xt:.3g}</text>')
    for i in range(5):
        yt = float(y_values[0] + i * (y_values[-1] - y_values[0]) / 4.0)
        py = margin_top + plot_height - i * plot_height / 4.0
        elements.append(f'<line x1="{margin_left - 6}" y1="{py:.2f}" x2="{margin_left}" y2="{py:.2f}" stroke="#222"/>')
        elements.append(f'<text x="{margin_left - 10}" y="{py + 4:.2f}" text-anchor="end" font-size="12" font-family="monospace">{yt:.3g}</text>')

    legend_x = width - 270
    legend_y = margin_top + 12
    elements.append(f'<rect x="{legend_x}" y="{legend_y}" width="220" height="{28 * len(labels) + 10}" fill="#fafafa" stroke="#bbb"/>')
    for idx, key in enumerate(sorted(labels)):
        y0 = legend_y + 20 + 26 * idx
        elements.append(f'<rect x="{legend_x + 12}" y="{y0 - 10}" width="18" height="18" fill="{palette[key]}" stroke="#444"/>')
        elements.append(f'<text x="{legend_x + 40}" y="{y0 + 3}" font-size="13" font-family="monospace">{labels[key]}</text>')

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
    ymax = float(np.max(values))
    if ymax <= 0.0:
        ymax = 1.0
    ymax *= 1.1
    bar_w = plot_width / max(1, len(values))

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width/2:.1f}" y="28" text-anchor="middle" font-size="20" font-family="monospace">{title}</text>',
        f'<line x1="{margin_left}" y1="{margin_top + plot_height}" x2="{margin_left + plot_width}" y2="{margin_top + plot_height}" stroke="#222" stroke-width="1.5"/>',
        f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_height}" stroke="#222" stroke-width="1.5"/>',
    ]
    for i in range(5):
        yt = i * ymax / 4.0
        py = margin_top + plot_height - i * plot_height / 4.0
        elements.append(f'<line x1="{margin_left - 6}" y1="{py:.2f}" x2="{margin_left}" y2="{py:.2f}" stroke="#222"/>')
        elements.append(f'<text x="{margin_left - 10}" y="{py + 4:.2f}" text-anchor="end" font-size="12" font-family="monospace">{yt:.3g}</text>')
        elements.append(f'<line x1="{margin_left}" y1="{py:.2f}" x2="{margin_left + plot_width}" y2="{py:.2f}" stroke="#ddd" stroke-width="1"/>')

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
