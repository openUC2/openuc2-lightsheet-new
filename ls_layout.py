"""
ls_layout.py - the light-sheet arrangements on the 50 mm cube grid: a drawing in the style of
the kit sketch (top view, one cube layer), the module of every cell with its CAD-new name, the
beam path, and the numbers from ls_optics.py. Writes layouts/<key>/ls_layout.{png,json}.

    uv run --with optiland --with numpy --with scipy --with matplotlib python ls_layout.py

Cells are (column, row), row 0 at the top as in the sketch; cell centres at (50 c, -50 r).
The sample hangs from the XYZ stage above the aquarium (not drawn).
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ls_optics as O  # noqa: E402

G = O.GRID
INK, INK2, WATER = "#0b0b0b", "#52514e", "#dbeafb"
BEAM, DET, SHEET = "#2a78d6", "#1baf7a", "#2a78d6"


@dataclass
class Cell:
    cell: Tuple[int, int]
    kind: str                 # fibre, lens+, lens-, cyl, mirror, kinmirror, obj4x, aquarium, obj10x, zdrive, filter, camera, free
    label: str
    module: str               # CAD-new / kit name, or what has to be made
    source: str               # kit it comes from, or "new"
    beam_dir: str = "+x"      # direction of the beam through the cell
    note: str = ""


@dataclass
class Arrangement:
    key: str
    title: str
    cols: int
    rows: int
    cells: List[Cell]
    detect_water_mm: float
    na: float = 0.25
    notes: List[str] = field(default_factory=list)


def arrangements() -> Dict[str, Arrangement]:
    det = [Cell((2, 1), "obj10x", "Electronic z-stage\nwith 10x imaging\nobjective", "L Plan 10x (InfinityBox) on the electronic z-stage", "LightSheet add-on", "-y"),
           Cell((1, 1), "zdrive", "", "z-stage drive (stepper, lead screw)", "LightSheet add-on", "-y"),
           Cell((1, 2), "zdrive2", "", "", "", "-y"),
           Cell((2, 2), "filter", "Emission filter", "ASS - 2053 - CUBREMOVAEMIFILWLS485 - V04 - A", "LightSheet add-on", "-y"),
           Cell((2, 3), "camera", "Camera with\ntube lens f100", "CLENS-100TEL + HIK MV-CE060-10UM-PRO (IMX178)", "InfinityBox", "-y")]
    k = Arrangement("K", "Today (kit, as documented)", 3, 4, [
        Cell((0, 3), "fibre", "Fibre + collimator\nf28", "fibre FC/PC (PRT-2161 INSFIBHOL) + f28 D16 lens (PRT-2102)", "FluoBox laser", "+y"),
        Cell((0, 2), "free", "(old f25 lens cube)", "-", "", "+y"),
        Cell((0, 1), "cyl", "Cylindrical lens\nf100 (power out of\nthe drawing)", "ASS - 2054 - CUBLENCYL30x32F100 - V04 - A", "LightSheet add-on", "+y", "z"),
        Cell((0, 0), "mirror", "Mirror", "ASS - 2020 - CUBMIR45°TH2 - V04", "CoreBox", "+y"),
        Cell((1, 0), "obj4x", "4x illuminating\nobjective", "4x objective in ASS - 2088 - CUBRMSTHRMOV - V04", "CoreBox", "+x"),
        Cell((2, 0), "aquarium", "Aquarium\n(sample chamber)", "ASS - 2052 - CUBAQU - V04", "LightSheet add-on", "+x"),
        Cell((1, 3), "free", "", "-", "", "+x")] + det, 9.08,
        notes=["Sheet 2.3 um thin but only 0.07 mm long; the camera sees 1.47 x 0.98 mm.",
               "Detection through 9.1 mm of water: Strehl 0.27 at NA 0.25."])
    a = Arrangement("A", "A - compact, same 3 x 4 footprint (recommended)", 3, 4, [
        Cell((0, 3), "fibre", "Fibre + collimator\nf28", "fibre FC/PC (PRT-2161 INSFIBHOL) + f28 D16 lens (PRT-2102)", "FluoBox laser", "+y"),
        Cell((0, 2), "lens+", "Lens f100", "ASS - 2022 - CUBLEND40F100 - V04", "CoreBox", "+y"),
        Cell((0, 1), "lens-", "Lens f-50\n(2:1 reducer)", "ASS - 2023 - CUBLEND43F-50 - V04", "CoreBox", "+y"),
        Cell((0, 0), "kinmirror", "Kinematic mirror", "ASS - 2013 - CUBKINMIR45°90° - V04", "CoreBox / kit", "+y"),
        Cell((1, 0), "cyl", "Cylindrical lens f50\n(power in the drawing)", "f50 cylinder in a rotatable cartridge", "new holder", "+x", "y"),
        Cell((2, 0), "aquarium", "Aquarium, detection\nwindow 3 mm from\nthe sheet", "new insert for ASS - 2052 - CUBAQU - V04", "new insert", "+x"),
        Cell((1, 3), "free", "free", "-", "", "+x")] + det, 3.0,
        notes=["4x objective no longer needed; CoreBox lenses f100 and f-50 shrink the beam to 2.4 mm.",
               "Cylinder 2.9 mm behind its cell centre (the water shifts the focus)."])
    a2 = Arrangement("A2", "A2 - fewest parts", 3, 4, [
        Cell((0, 1), "fibre", "Fibre collimator\nf11.2", "fibre collimator f11.2 (BUY DM11x27.8) in a round clamp", "new holder", "+y"),
        Cell((0, 2), "free", "free", "-", "", "+y"),
        Cell((0, 3), "free", "free", "-", "", "+y"),
        Cell((0, 0), "kinmirror", "Kinematic mirror", "ASS - 2013 - CUBKINMIR45°90° - V04", "CoreBox / kit", "+y"),
        Cell((1, 0), "cyl", "Cylindrical lens f50\n(power in the drawing)", "f50 cylinder in a rotatable cartridge", "new holder", "+x", "y"),
        Cell((2, 0), "aquarium", "Aquarium, detection\nwindow 3 mm from\nthe sheet", "new insert for ASS - 2052 - CUBAQU - V04", "new insert", "+x"),
        Cell((1, 3), "free", "free", "-", "", "+x")] + det, 3.0,
        notes=["Short collimator gives a 1.9 mm beam directly; sheet 1.9 mm tall (crop the long side)."])
    det_b = [replace_cell(c, (c.cell[0] + 1, c.cell[1])) for c in det]
    b = Arrangement("B", "B - tall, uniform sheet, 4 x 4", 4, 4, [
        Cell((0, 1), "fibre", "Fibre + collimator\nf28", "fibre FC/PC (PRT-2161 INSFIBHOL) + f28 D16 lens (PRT-2102)", "FluoBox laser", "+y"),
        Cell((0, 0), "kinmirror", "Kinematic mirror", "ASS - 2013 - CUBKINMIR45°90° - V04", "CoreBox / kit", "+y"),
        Cell((1, 0), "cyl", "Cylindrical lens f100\n(turned: power in\nthe drawing)", "ASS - 2054 - CUBLENCYL30x32F100 - V04 - A, turned 90 deg", "LightSheet add-on", "+x", "y"),
        Cell((2, 0), "free", "free (later: slit\nor pivot scanner)", "-", "", "+x"),
        Cell((3, 0), "aquarium", "Aquarium, detection\nwindow 3 mm from\nthe sheet", "new insert for ASS - 2052 - CUBAQU - V04", "new insert", "+x"),
        Cell((0, 2), "free", "free", "-", "", "+x"), Cell((0, 3), "free", "free", "-", "", "+x"),
        Cell((1, 1), "free", "free", "-", "", "+x"), Cell((1, 2), "free", "free", "-", "", "+x"), Cell((1, 3), "free", "free", "-", "", "+x"),
        Cell((2, 3), "free", "free", "-", "", "+x")] + det_b, 3.0,
        notes=["Kit's f100 cylinder turned 90 deg, two cells before the sheet; no extra lenses.",
               "Sheet 4.8 mm tall: 83 % brightness at the long edge of the field."])
    return {"K": k, "A": a, "A2": a2, "B": b}


def replace_cell(c: Cell, cell: Tuple[int, int]) -> Cell:
    return Cell(cell, c.kind, c.label, c.module, c.source, c.beam_dir, c.note)


def centre(cell: Tuple[int, int]) -> Tuple[float, float]:
    return (G * cell[0], -G * cell[1])


# ----------------------------------------------------------------------------- pictograms
def _rot(pts, d: str):
    """Pictograms are drawn for a beam along +x; turn them for the cell's beam direction."""
    ang = {"+x": 0, "+y": 90, "-x": 180, "-y": 270}[d]
    a = math.radians(ang)
    R = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    return np.asarray(pts, float) @ R.T


def _poly(ax, pts, d, cx, cy, **kw):
    from matplotlib.patches import Polygon
    p = _rot(pts, d) + (cx, cy)
    ax.add_patch(Polygon(p, closed=kw.pop("closed", True), **kw))


def _line(ax, pts, d, cx, cy, **kw):
    p = _rot(pts, d) + (cx, cy)
    ax.plot(p[:, 0], p[:, 1], **kw)


def draw_cell(ax, c: Cell):
    cx, cy = centre(c.cell)
    d = c.beam_dir
    lw = 1.6
    k = c.kind
    if k == "lens+":
        t = np.linspace(-1, 1, 30)
        pts = [(3.2 * (1 - tt ** 2) ** 0.5 * 1.0, 16 * tt) for tt in t] + [(-3.2 * (1 - tt ** 2) ** 0.5, 16 * tt) for tt in t[::-1]]
        _poly(ax, pts, d, cx, cy, facecolor="white", edgecolor=INK, lw=lw)
    elif k == "lens-":
        t = np.linspace(-1, 1, 30)
        pts = [(1.0 + 2.6 * tt ** 2, 16 * tt) for tt in t] + [(-1.0 - 2.6 * tt ** 2, 16 * tt) for tt in t[::-1]]
        _poly(ax, pts, d, cx, cy, facecolor="white", edgecolor=INK, lw=lw)
    elif k == "cyl":
        _poly(ax, [(-3, -15), (3, -15), (3, 15), (-3, 15)], d, cx, cy, facecolor="white", edgecolor=INK, lw=lw)
        for y in np.linspace(-12, 12, 9):
            _line(ax, [(-3, y), (3, y)] if c.note == "z" else [(-2.2, y), (2.2, y)], d, cx, cy, color=INK, lw=0.9)
        if c.note == "y":            # power in the drawing plane: show the curved face towards the source
            t = np.linspace(-1, 1, 30)
            _line(ax, [(-3 - 2.5 * (1 - tt ** 2), 15 * tt) for tt in t], d, cx, cy, color=INK, lw=lw)
    elif k in ("mirror", "kinmirror"):
        # 45 deg mirror folding +y -> +x drawn in the cell frame directly
        ax.plot([cx - 16, cx + 16], [cy - 16, cy + 16], color=INK, lw=2.4)
        for s in np.linspace(-14, 14, 8):
            ax.plot([cx + s, cx + s - 3], [cy + s, cy + s + 3], color=INK, lw=0.9)
        if k == "kinmirror":
            for s in (-9, 9):
                ax.add_patch(__import__("matplotlib.patches", fromlist=["Circle"]).Circle((cx + s - 9, cy + s + 9), 2.2, facecolor="white", edgecolor=INK, lw=1.2))
    elif k == "obj4x":
        _poly(ax, [(-14, -11), (6, -11), (6, 11), (-14, 11)], d, cx, cy, facecolor="white", edgecolor=INK, lw=lw)
        t = np.linspace(-math.pi / 2, math.pi / 2, 30)
        _line(ax, [(6 + 7 * math.cos(tt), 7 * math.sin(tt)) for tt in t], d, cx, cy, color=INK, lw=lw)
    elif k == "obj10x":
        _poly(ax, [(-9, -11), (9, -11), (9, 11), (-9, 11)], "+x", cx, cy, facecolor="white", edgecolor=INK, lw=lw)
        t = np.linspace(0, math.pi, 30)
        _line(ax, [(8 * math.cos(tt), 11 + 6 * math.sin(tt)) for tt in t], "+x", cx, cy, color=INK, lw=lw)
        ax.add_patch(__import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle((cx - 9 - 30, cy - 14), 30 + 9, 5, facecolor=INK, edgecolor=INK))
    elif k == "zdrive":
        ax.add_patch(__import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle((cx - 9, cy - 4), 18, 18, facecolor=INK, edgecolor=INK))
        ax.plot([cx, cx], [cy - 4, cy - 25], color=INK, lw=2)
    elif k == "zdrive2":
        ys = np.linspace(cy + 25, cy - 12, 10)
        ax.plot([cx, cx], [cy + 25, cy - 14], color=INK, lw=2)
        for y in ys[1:-1]:
            ax.plot([cx - 5, cx + 5], [y - 2, y + 2], color=INK, lw=1.4)
    elif k == "filter":
        _line(ax, [(0, -15), (0, -3)], "+y", cx, cy, color=INK, lw=3)
        _line(ax, [(0, 3), (0, 15)], "+y", cx, cy, color=INK, lw=3)
    elif k == "camera":
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((cx - 12, cy - 20), 24, 22, facecolor=INK, edgecolor=INK))
        ax.add_patch(Rectangle((cx - 7, cy + 2), 14, 10, facecolor="white", edgecolor=INK, lw=1.4))
        t = np.linspace(-1, 1, 20)
        ax.plot(cx + 8 * t, cy + 14 + 2.4 * np.sqrt(1 - t ** 2), color=INK, lw=1.2)
        ax.plot(cx + 8 * t, cy + 14 - 2.4 * np.sqrt(1 - t ** 2), color=INK, lw=1.2)
    elif k == "fibre":
        from matplotlib.patches import Rectangle
        _poly(ax, [(-20, -5), (-10, -5), (-10, 5), (-20, 5)], d, cx, cy, facecolor=INK, edgecolor=INK)
        t = np.linspace(-1, 1, 20)
        _poly(ax, [(8 + 0 * tt, 8 * tt) for tt in t] + [(8 + 2.2 * (1 - tt ** 2) ** 0.5, 8 * tt) for tt in t[::-1]], d, cx, cy,
              facecolor="white", edgecolor=INK, lw=1.4)
        tail = _rot([(-20, 0), (-24, -3), (-28, 3), (-32, 0)], d) + (cx, cy)
        ax.plot(tail[:, 0], tail[:, 1], color="#eda100", lw=2)


def draw_aquarium(ax, cell, d_window: float, sheet_y: float = 0.0):
    from matplotlib.patches import Circle, Rectangle
    cx, cy = centre(cell)
    ax.add_patch(Rectangle((cx - 17, cy - d_window), 34, 17 + d_window, facecolor=WATER, edgecolor="none"))
    ax.plot([cx - 17, cx + 17], [cy - d_window, cy - d_window], color=INK2, lw=2.2)          # detection window
    ax.plot([cx - 11.08, cx - 11.08], [cy - d_window, cy + 17], color=INK2, lw=2.2)           # illumination window
    ax.add_patch(Circle((cx + 0.0, cy + sheet_y), 2.35, facecolor="none", edgecolor=INK, lw=1.0, ls="--"))
    ax.text(cx + 3, cy + 7, "sample", fontsize=6.5, color=INK2)


def draw_beams(ax, arr: Arrangement, sh: "O.Sheet"):
    """Illumination as a band (1/e2 width x 3 for visibility), the sheet focus, the detection axis."""
    cells = {c.kind: c for c in arr.cells}
    by_cell = {c.cell: c for c in arr.cells}
    fib = next(c for c in arr.cells if c.kind == "fibre")
    mir = next(c for c in arr.cells if c.kind in ("mirror", "kinmirror"))
    aq = next(c for c in arr.cells if c.kind == "aquarium")
    foc = next(c for c in arr.cells if c.kind in ("cyl", "obj4x") and c.beam_dir == "+x") if arr.key == "K" else next(c for c in arr.cells if c.kind == "cyl")
    x0, y0 = centre(fib.cell)
    xm, ym = centre(mir.cell)
    xa, ya = centre(aq.cell)
    xf, yf = centre(foc.cell)
    scale = 3.0
    v = O.variants()[arr.key]
    # beam half-widths along the vertical leg (in the drawing plane = thickness axis)
    def hw(s):
        return scale * O.w_of(O.q_at(v, "y", s, sh.stations))
    s_m = -(xa - xm)          # path position of the mirror (from the sheet)
    pts_up = []
    for y in np.linspace(y0 + 8, ym, 40):
        s = s_m - (ym - y)
        pts_up.append((y, hw(s)))
    left = [(x0 - w, y) for y, w in pts_up]
    right = [(x0 + w, y) for y, w in pts_up]
    from matplotlib.patches import Polygon
    ax.add_patch(Polygon(left + right[::-1], closed=True, facecolor=BEAM, alpha=0.35, edgecolor="none"))
    pts_x = []
    for x in np.linspace(xm, xa, 120):
        s = -(xa - x)
        w = O.w_of(O.q_at(v, "y", s, sh.stations))
        pts_x.append((x, max(scale * w, 0.25)))
    top = [(x, ya + w) for x, w in pts_x]
    bot = [(x, ya - w) for x, w in pts_x]
    ax.add_patch(Polygon(top + bot[::-1], closed=True, facecolor=BEAM, alpha=0.35, edgecolor="none"))
    ax.plot([xa - 17, xa + 17], [ya, ya], color=SHEET, lw=1.0)
    # detection axis
    cam = next(c for c in arr.cells if c.kind == "camera")
    xc, yc = centre(cam.cell)
    ax.plot([xa, xa], [ya, yc + 12], color=DET, lw=2.0, ls=(0, (5, 3)))


def draw(arr: Arrangement, out: Path, numbers: Dict) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    v = O.variants()[arr.key]
    sh = O.sheet(v)
    fig_w = 4.2 + 2.2 * arr.cols
    fig, ax = plt.subplots(figsize=(fig_w, 8.6))
    ax.set_aspect("equal")
    ax.axis("off")
    for c in range(arr.cols):
        for r in range(arr.rows):
            cx, cy = centre((c, r))
            ax.add_patch(Rectangle((cx - G / 2, cy - G / 2), G, G, facecolor="white", edgecolor=INK, lw=2.2))
    aq = next(c for c in arr.cells if c.kind == "aquarium")
    draw_aquarium(ax, aq.cell, arr.detect_water_mm)
    draw_beams(ax, arr, sh)
    for c in arr.cells:
        draw_cell(ax, c)
    # labels: column 0 to the left, last column to the right, middle columns above (row 0) or inside
    for c in arr.cells:
        if not c.label:
            continue
        cx, cy = centre(c.cell)
        col, row = c.cell
        if c.kind == "free":
            ax.text(cx, cy - (15 if row == 0 else 0), c.label, ha="center", va="center", fontsize=7, color="#9a9890")
        elif col == 0:
            ax.text(cx - G / 2 - 4, cy, c.label, ha="right", va="center", fontsize=8.5, color=INK)
        elif col == arr.cols - 1:
            ax.text(cx + G / 2 + 4, cy, c.label, ha="left", va="center", fontsize=8.5, color=INK)
        elif row == 0:
            ax.text(cx, cy + G / 2 + 4, c.label, ha="center", va="bottom", fontsize=8.5, color=INK)
        else:
            ax.text(cx, cy - G / 2 + 3, c.label, ha="center", va="bottom", fontsize=7, color=INK)
    # numbers box
    q = numbers["detection"]
    nb = [f"sheet: waist FWHM {sh.fwhm0_um:.1f} um, length {2 * sh.zr_mm:.2f} mm, height {sh.height_mm:.1f} mm",
          f"sheet FWHM at the field's short edge (0.49 mm): {float(sh.fwhm_at(0.49)):.1f} um",
          f"detection: {arr.detect_water_mm:.1f} mm water, NA {arr.na:.2f}: Strehl {q['strehl']:.2f}, "
          f"2 um modulation {q['mtf_2um']:.2f} (perfect 0.37)",
          "camera: 5x, 0.48 um pixels, field 1.47 x 0.98 mm",
          "beam width drawn 3x; the sheet waist sits on the aquarium centre"] + arr.notes
    x_left = -G / 2
    ax.text(x_left, -G * (arr.rows - 0.5) - 16, "\n".join(nb), ha="left", va="top", fontsize=7.6, color=INK2, linespacing=1.45)
    ax.set_title(arr.title, fontsize=11, color=INK, loc="left")
    xmin = -G / 2 - 95
    xmax = G * (arr.cols - 0.5) + 95
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-G * (arr.rows - 0.5) - 96, G / 2 + 26)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def layout_json(arr: Arrangement, numbers: Dict) -> Dict:
    v = O.variants()[arr.key]
    sh = O.sheet(v)
    return {
        "key": arr.key, "title": arr.title, "grid_mm": G, "cols": arr.cols, "rows": arr.rows,
        "frame": "cells (column, row), row 0 at the top; cell centre (50 c, -50 r); one cube layer; sample from above on the XYZ stage",
        "cells": [{"cell": list(c.cell), "kind": c.kind, "label": c.label.replace("\n", " "), "module": c.module, "source": c.source,
                   "beam_direction": c.beam_dir} for c in arr.cells if c.kind != "zdrive2"],
        "lenses_on_path": [{"lens": st.lens.name, "part": st.lens.part, "cell": list(st.cell), "path_from_sheet_mm": st.s,
                            "moved_from_cell_centre_mm": st.s - st0.s} for st, st0 in zip(sh.stations, v.stations)],
        "sheet": {"waist_fwhm_um": sh.fwhm0_um, "length_2zr_mm": 2 * sh.zr_mm, "height_mm": sh.height_mm,
                  "fwhm_at_0.49_um": float(sh.fwhm_at(0.49)), "fwhm_at_0.74_um": float(sh.fwhm_at(0.74)),
                  "irradiance_w_cm2_per_mw": sh.peak_irradiance_w_cm2_per_mw()},
        "detection": dict(numbers["detection"], detect_water_mm=arr.detect_water_mm),
        "notes": arr.notes,
    }


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    out_root = HERE / "layouts"
    arrs = arrangements()
    paths = []
    for key, arr in arrs.items():
        dq = O.detection_quality(arr.detect_water_mm, O.Detection(na=arr.na))
        numbers = {"detection": {"strehl": dq.strehl, "mtf_2um": dq.mtf_2um, "mtf_4um": dq.mtf_4um, "rms_waves": dq.rms_waves,
                                 "axial_fwhm_um": dq.axial_fwhm_um, "na": arr.na}}
        p = draw(arr, out_root / key / "ls_layout.png", numbers)
        (out_root / key / "ls_layout.json").write_text(json.dumps(layout_json(arr, numbers), indent=2, ensure_ascii=False), encoding="utf-8")
        paths.append(p)
        print("wrote", p)
    # overview: the four drawings side by side
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    imgs = [plt.imread(str(p)) for p in paths]
    fig, axs = plt.subplots(1, len(imgs), figsize=(6 * len(imgs), 8.5), gridspec_kw={"width_ratios": [im.shape[1] / im.shape[0] for im in imgs]})
    for ax, im in zip(axs, imgs):
        ax.imshow(im)
        ax.axis("off")
    fig.tight_layout()
    ov = HERE / "docs" / "ls_arrangements.png"
    fig.savefig(ov, dpi=110, facecolor="white", bbox_inches="tight")
    print("wrote", ov)


if __name__ == "__main__":
    main()
