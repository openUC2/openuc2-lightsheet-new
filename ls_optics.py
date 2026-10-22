"""
ls_optics.py - optics of the openUC2 light-sheet microscope.

The sheet is a Gaussian beam traced per axis (ABCD): y = sheet thickness = detection axis,
z = sheet height, x = propagation. Optiland models of the illumination (unfolded) and of
the detection through the aquarium window check the Gaussian numbers and give the
aberrations of imaging into water with an air objective.

    uv run --with optiland --with numpy --with scipy --with matplotlib python ls_optics.py --validate
    ... python ls_optics.py --report docs/ls_numbers.md --figures docs --export layouts

Units: mm for lengths, um for wavelengths and spot sizes unless a name says otherwise.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent

LAM_EX = 0.488          # um, laser
LAM_EM = 0.520          # um, emission (GFP, fluorescein)
N_WATER_EX = 1.3371     # water at 488 nm, 20 C
N_WATER_EM = 1.3354     # water at 520 nm
N_COVER = 1.523         # cover glass (D 263 M class) at 500 nm
N_BK7_EX = 1.5224       # N-BK7 at 488 nm
GRID = 50.0             # cube pitch (D-E1)


# ----------------------------------------------------------------------------- parts
@dataclass(frozen=True)
class Fibre:
    """Single-mode fibre, 4 um core; the mode-field diameter at 488 nm is assumed."""
    mfd_um: float = 3.6

    @property
    def w_mm(self) -> float:
        return self.mfd_um / 2e3

    def theta(self, lam_um: float = LAM_EX) -> float:
        """1/e^2 half-angle of the far field (rad)."""
        return lam_um * 1e-3 / (math.pi * self.w_mm)


@dataclass(frozen=True)
class Lens:
    """Lens on the unfolded path. axes: 'yz' spherical, 'y' cylinder with power across the
    sheet (thickness), 'z' cylinder with power along the sheet height. R/ct given = real
    plano-convex (curved side to the source) in the Optiland model, else paraxial."""
    name: str
    f: float
    axes: str = "yz"
    part: str = ""
    R: Optional[float] = None
    ct: float = 0.0
    glass: str = "N-BK7"

    @property
    def n(self) -> float:
        return N_BK7_EX


def cyl_lens(name: str, f: float, axes: str, part: str, ct: float = 5.0) -> Lens:
    """Plano-convex N-BK7 cylinder with the radius that gives f at 488 nm."""
    return Lens(name, f, axes, part, R=f * (N_BK7_EX - 1.0), ct=ct)


COLLIMATOR_F28 = Lens("collimator f28 D16 (FluoBox laser)", 28.0, part="BUY - Lens - f28 D16 sI3.6 R20 pl-cx + PRT-2102 MASINSLEND16F28")
COLLIMATOR_F11 = Lens("fibre collimator f11.2", 11.2, part="BUY fibre collimator DM11x27.8 f11.2")
LENS_F100 = Lens("lens f100 D40 (CoreBox)", 100.0, part="ASS - 2022 - CUBLEND40F100 - V04")
LENS_FM50 = Lens("lens f-50 D43 (CoreBox)", -50.0, part="ASS - 2023 - CUBLEND43F-50 - V04")
CYL_F100_Z = cyl_lens("cylinder f100, power along the sheet height", 100.0, "z", "ASS - 2054 - CUBLENCYL30x32F100 - V04 - A")
CYL_F100_Y = cyl_lens("cylinder f100, power across the sheet", 100.0, "y", "ASS - 2054 - CUBLENCYL30x32F100 - V04 - A (turned 90 deg)")
CYL_F50_Y = cyl_lens("cylinder f50, power across the sheet", 50.0, "y", "user's f50 cylinder (holder to be made)")


def objective_4x(f: float) -> Lens:
    return Lens(f"4x objective (f {f:g})", f, part="BUY - Objectiv 4x - XXX in ASS - 2088 - CUBRMSTHRMOV - V04")


@dataclass(frozen=True)
class Chamber:
    """ASS-2052 CUBAQU: 0.17 mm cover-glass windows; sheet on the cube centre."""
    t_glass: float = 0.17
    illum_water: float = 11.08     # illumination window (inner face) -> sheet
    detect_water: float = 9.08     # detection window (inner face) -> sheet


@dataclass(frozen=True)
class Detection:
    """L Plan 10x (infinity) on the z-stage, CLENS-100TEL f100 tube lens, HIK MV-CE060-10UM-PRO."""
    na: float = 0.25
    f_obj: float = 20.0            # 200 mm design tube length assumed
    f_tl: float = 100.0
    wd_air: float = 16.3           # parfocal 45 - nose 28.7 (CAD), to be confirmed
    pixel_um: float = 2.4
    nx: int = 3072
    ny: int = 2048
    obj_to_tl: float = 100.0       # infinity space (objective cell -> camera cell)
    tl_clear_r: float = 11.0       # clear radius of filter / tube lens

    @property
    def mag(self) -> float:
        return self.f_tl / self.f_obj

    @property
    def px_obj_um(self) -> float:
        return self.pixel_um / self.mag

    @property
    def fov_mm(self) -> Tuple[float, float]:
        return (self.nx * self.px_obj_um * 1e-3, self.ny * self.px_obj_um * 1e-3)


# ----------------------------------------------------------------------------- variants
@dataclass(frozen=True)
class Station:
    s: float                      # path position of the lens (thin lens / rear principal plane), from the sheet
    lens: Lens
    cell: Tuple[int, int]
    solve: bool = False           # moved along the path so the sheet waist sits on the sheet plane


@dataclass(frozen=True)
class Variant:
    key: str
    title: str
    collimator: Lens
    s_collimator: float
    stations: Tuple[Station, ...]
    fibre_cell: Tuple[int, int]
    chamber: Chamber = Chamber()
    fibre: Fibre = Fibre()
    note: str = ""


def variants(f_4x: float = 31.0) -> Dict[str, Variant]:
    """The arrangements compared in round 1 (cells as (column, row), row 0 = top)."""
    o4 = objective_4x(f_4x)
    return {
        "K": Variant("K", f"kit as documented: f28 collimator, f100 cylinder into the pupil of the 4x (f {f_4x:g})",
                     COLLIMATOR_F28, -250.0,
                     (Station(-150.0, CYL_F100_Z, (0, 1)), Station(-50.0, o4, (1, 0), solve=True)), (0, 3),
                     note="cylinder placed with its line focus in the 4x pupil, as in the alignment guide"),
        "A": Variant("A", "compact 3 x 4: f28 collimator, 2:1 Galilean reducer (f100, f-50), f50 cylinder",
                     COLLIMATOR_F28, -250.0,
                     (Station(-200.0, LENS_F100, (0, 2)), Station(-150.0, LENS_FM50, (0, 1)),
                      Station(-50.0, CYL_F50_Y, (1, 0), solve=True)), (0, 3)),
        "A0": Variant("A0", "A without the reducer: thin, short sheet for small regions",
                      COLLIMATOR_F28, -250.0, (Station(-50.0, CYL_F50_Y, (1, 0), solve=True),), (0, 3)),
        "A2": Variant("A2", "fewest parts: f11.2 fibre collimator, f50 cylinder",
                      COLLIMATOR_F11, -150.0, (Station(-50.0, CYL_F50_Y, (1, 0), solve=True),), (0, 1)),
        "B": Variant("B", "4 x 4, tall uniform sheet: f28 collimator, f100 cylinder two cells before the sheet",
                     COLLIMATOR_F28, -200.0, (Station(-100.0, CYL_F100_Y, (1, 0), solve=True),), (0, 1)),
    }


# ----------------------------------------------------------------------------- Gaussian sheet (ABCD)
def reduced(s1: float, s2: float, ch: Chamber, n_w: float = N_WATER_EX) -> float:
    """Reduced length (sum of t/n) from path position s1 to s2 (air, window glass, water)."""
    if s2 <= s1:
        return 0.0
    edges = [(-math.inf, -ch.illum_water - ch.t_glass, 1.0), (-ch.illum_water - ch.t_glass, -ch.illum_water, N_COVER),
             (-ch.illum_water, math.inf, n_w)]
    total = 0.0
    for a, b, n in edges:
        lo, hi = max(s1, a), min(s2, b)
        if hi > lo:
            total += (hi - lo) / n
    return total


def _events(v: Variant, stations: Sequence[Station]) -> List[Tuple[float, Lens]]:
    return sorted([(v.s_collimator, v.collimator)] + [(st.s, st.lens) for st in stations], key=lambda e: e[0])


def q_at(v: Variant, axis: str, s_end: float, stations: Optional[Sequence[Station]] = None, lam: float = LAM_EX) -> complex:
    """Reduced complex beam parameter at path position s_end for one axis ('y' or 'z')."""
    stations = v.stations if stations is None else stations
    lam_mm = lam * 1e-3
    q = 1j * math.pi * v.fibre.w_mm ** 2 / lam_mm
    s = v.s_collimator - v.collimator.f
    for s_e, lens in _events(v, stations):
        if s_e > s_end:
            break
        q += reduced(s, s_e, v.chamber)
        s = s_e
        if axis in lens.axes:
            q = 1.0 / (1.0 / q - 1.0 / lens.f)
    return q + reduced(s, s_end, v.chamber)


def w_of(q: complex, lam: float = LAM_EX) -> float:
    """1/e^2 radius (mm) from a reduced q."""
    return math.sqrt(-lam * 1e-3 / (math.pi * (1.0 / q).imag))


def solve_stations(v: Variant) -> Tuple[Station, ...]:
    """Move the 'solve' station so the y-waist lands on the sheet (s = 0); K's cylinder follows
    its objective (line focus in the pupil)."""
    sts = list(v.stations)
    for i, st in enumerate(sts):
        if not st.solve:
            continue

        def waist_offset(s_try: float) -> float:
            trial = sts[:i] + [replace(st, s=s_try)] + sts[i + 1:]
            q_i = q_at(v, "y", -v.chamber.illum_water, trial)
            return -N_WATER_EX * q_i.real - v.chamber.illum_water
        lo, hi = st.s - 30.0, st.s + 30.0
        f_lo, f_hi = waist_offset(lo), waist_offset(hi)
        if f_lo * f_hi > 0:
            raise ValueError(f"{v.key}: cannot put the waist on the sheet by moving {st.lens.name}")
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f_mid = waist_offset(mid)
            if f_lo * f_mid <= 0:
                hi = mid
            else:
                lo, f_lo = mid, f_mid
        sts[i] = replace(st, s=0.5 * (lo + hi))
    if v.key.startswith("K"):
        obj = next(st for st in sts if st.solve)
        cyl = next(j for j, st in enumerate(sts) if st.lens.axes == "z")
        sts[cyl] = replace(sts[cyl], s=obj.s - obj.lens.f - sts[cyl].lens.f)
    return tuple(sts)


@dataclass
class Sheet:
    key: str
    beam_d_mm: float              # 1/e^2 diameter of the collimated beam at the focusing lens
    na_air: float                 # 1/e^2 half-angle of the thickness axis, in air
    w0_um: float
    fwhm0_um: float
    zr_mm: float                  # Rayleigh range in water
    waist_offset_mm: float        # waist - sheet plane (should be 0 after solving)
    height_mm: float              # 1/e^2 full height at the sheet
    stations: Tuple[Station, ...]
    lens_shift_mm: Dict[str, float]

    def w_at(self, x_mm) -> np.ndarray:
        x = np.asarray(x_mm, float)
        return self.w0_um * np.sqrt(1.0 + ((x - self.waist_offset_mm) / self.zr_mm) ** 2)

    def fwhm_at(self, x_mm) -> np.ndarray:
        return self.w_at(x_mm) * math.sqrt(2.0 * math.log(2.0))

    def height_fraction(self, z_mm: float) -> float:
        return math.exp(-2.0 * (z_mm / (self.height_mm / 2.0)) ** 2)

    def peak_irradiance_w_cm2_per_mw(self) -> float:
        return 2.0 * 1e-3 / (math.pi * (self.w0_um * 1e-4) * (self.height_mm / 20.0))


def sheet(v: Variant) -> Sheet:
    sts = solve_stations(v)
    ch = v.chamber
    q_i = q_at(v, "y", -ch.illum_water, sts)
    x_w = -N_WATER_EX * q_i.real
    w0 = math.sqrt(LAM_EX * 1e-3 * q_i.imag / math.pi)
    zr = N_WATER_EX * q_i.imag
    q_z = q_at(v, "z", 0.0, sts)
    focus = next(st for st in sts if st.solve)
    s_before = max(focus.s - 1e-6, v.s_collimator + 1e-6)
    beam_d = 2.0 * w_of(q_at(v, "y", s_before, sts))
    shifts = {st.lens.name: st.s - st0.s for st, st0 in zip(sts, v.stations) if abs(st.s - st0.s) > 1e-9}
    return Sheet(v.key, beam_d, beam_d / 2.0 / abs(focus.lens.f), w0 * 1e3, w0 * 1e3 * math.sqrt(2 * math.log(2)), zr,
                 x_w - ch.illum_water, 2.0 * w_of(q_z), sts, shifts)


# ----------------------------------------------------------------------------- Optiland: illumination
def _material(m):
    from optiland.materials import IdealMaterial
    if isinstance(m, str) and m.startswith("ideal:"):
        return IdealMaterial(n=float(m.split(":")[1]))
    return m


def illumination_specs(v: Variant, stations: Optional[Sequence[Station]] = None) -> List[Tuple[int, dict]]:
    """Optiland surfaces of the unfolded illumination path (fibre -> sheet plane)."""
    sts = solve_stations(v) if stations is None else stations
    rows: List[Tuple[float, dict]] = []
    for s, lens in _events(v, sts):
        if lens.R is None:
            rows.append((s, dict(surface_type="paraxial", f=lens.f, comment=lens.name)))
            continue
        s_front = s - lens.ct * (lens.n - 1.0) / lens.n          # rear principal plane on s
        if lens.axes == "y":
            kw = dict(surface_type="toroidal", radius_x=math.inf, radius_y=lens.R)
        elif lens.axes == "z":
            kw = dict(surface_type="toroidal", radius_x=lens.R, radius_y=math.inf)
        else:
            kw = dict(radius=lens.R)
        rows.append((s_front, dict(**kw, material=lens.glass, comment=f"{lens.name}, curved face")))
        rows.append((s_front + lens.ct, dict(radius=math.inf, comment=f"{lens.name}, flat face")))
    ch = v.chamber
    rows.append((-ch.illum_water - ch.t_glass, dict(radius=math.inf, material=f"ideal:{N_COVER}", comment="window, outer face")))
    rows.append((-ch.illum_water, dict(radius=math.inf, material=f"ideal:{N_WATER_EX}", comment="window, inner face (water)")))
    rows.append((0.0, dict(radius=math.inf, comment="sheet plane (aquarium cube centre)")))
    specs = [(0, dict(radius=math.inf, thickness=v.collimator.f, comment="fibre tip"))]
    for i, (s, kw) in enumerate(rows):
        kw = dict(kw)
        if i + 1 < len(rows):
            kw["thickness"] = rows[i + 1][0] - s
        if i == 0:
            kw["is_stop"] = True
        specs.append((i + 1, kw))
    return specs


def build(specs, aperture: Tuple[str, float], lam: float, fields: Sequence[float] = (0.0,), field_type: str = "angle"):
    from optiland import optic
    o = optic.Optic()
    for idx, kw in specs:
        kw = {k: (_material(v) if k == "material" else v) for k, v in kw.items()}
        o.surfaces.add(index=idx, **kw)
    o.set_aperture(aperture_type=aperture[0], value=aperture[1])
    o.set_field_type(field_type)
    for f in fields:
        o.add_field(y=f)
    o.add_wavelength(unit="um", value=lam, is_primary=True)
    return o


def illumination_optic(v: Variant):
    return build(illumination_specs(v), ("objectNA", math.sin(v.fibre.theta())), LAM_EX)


def trace_illumination(v: Variant, n_rays: int = 11) -> Dict[str, float]:
    """Real rays: focus of the thickness fan (in water, from the sheet plane), its rms spot,
    and the half-height of the height fan at the sheet plane."""
    o = illumination_optic(v)
    p = np.linspace(-1.0, 1.0, n_rays)
    o.trace_generic(0.0, 0.0, np.zeros(n_rays), p, LAM_EX)
    sg = o.surfaces
    y, M, N = np.asarray(sg.y)[-1], np.asarray(sg.M)[-1], np.asarray(sg.N)[-1]
    t = M / N
    z_f = -float(np.sum(y * t) / np.sum(t * t))
    spot = float(np.sqrt(np.mean((y + z_f * t) ** 2))) * 1e3
    o.trace_generic(0.0, 0.0, p, np.zeros(n_rays), LAM_EX)
    x = np.asarray(o.surfaces.x)[-1]
    return {"focus_mm": z_f, "rms_spot_um": spot, "half_height_mm": float(np.max(np.abs(x)))}


# ----------------------------------------------------------------------------- Optiland: detection
FRINGE_RADIAL = {   # rotationally symmetric fringe Zernikes: index -> (n, polynomial coefficients high->low in rho)
    4: (2, [2, 0, -1]),
    9: (4, [6, 0, -6, 0, 1]),
    16: (6, [20, 0, -30, 0, 12, 0, -1]),
    25: (8, [70, 0, -140, 0, 90, 0, -20, 0, 1]),
    36: (10, [252, 0, -630, 0, 560, 0, -210, 0, 30, 0, -1]),
    37: (12, [924, 0, -2772, 0, 3150, 0, -1680, 0, 420, 0, -42, 0, 1]),
}


def detection_specs(d_water: float, det: Detection, t_glass: float = 0.17) -> List[Tuple[int, dict]]:
    """Object point in water, window, air gap to a perfect objective (paraxial), tube lens, camera."""
    d = max(d_water, 1e-3)
    gap = det.f_obj - t_glass / N_COVER - d / N_WATER_EM
    return [
        (0, dict(radius=math.inf, thickness=d, material=f"ideal:{N_WATER_EM}", comment="sample point in water")),
        (1, dict(radius=math.inf, thickness=t_glass, material=f"ideal:{N_COVER}", comment="detection window")),
        (2, dict(radius=math.inf, thickness=gap, comment="air gap to the objective")),
        # 60 mm instead of the real 100 mm: with the pupil in the tube lens's front focus the image side is
        # telecentric and Optiland's chief-ray wavefront fails; on axis the gap does not change the wavefront
        (3, dict(surface_type="paraxial", f=det.f_obj, thickness=60.0, is_stop=True, comment="10x objective (perfect)")),
        (4, dict(surface_type="paraxial", f=det.f_tl, thickness=det.f_tl, comment="tube lens f100")),
        (5, dict(radius=math.inf, comment="camera")),
    ]


def detection_optic(d_water: float, det: Detection, t_glass: float = 0.17, field_mm: float = 0.0):
    """The aperture is the objective's pupil (stop size iterated to the NA in air): Optiland's
    chief-ray wavefront does not work with an object-NA aperture here."""
    fields = (0.0,) if field_mm == 0.0 else (0.0, field_mm)
    stop_d = 2.0 * det.f_obj * math.tan(math.asin(det.na))
    for _ in range(3):
        o = build(detection_specs(d_water, det, t_glass), ("float_by_stop_size", stop_d), LAM_EM, fields, "object_height")
        o.trace_generic(0.0, 0.0, np.array([0.0]), np.array([1.0]), LAM_EM)
        M, N = float(np.asarray(o.surfaces.M)[2][0]), float(np.asarray(o.surfaces.N)[2][0])
        na_air = abs(M) / math.hypot(M, N)
        stop_d *= det.na / na_air
    return o


def zernike_symmetric(d_water: float, det: Detection, t_glass: float = 0.17) -> Dict[int, float]:
    """Rotationally symmetric fringe Zernike coefficients (waves) of the on-axis wavefront."""
    from optiland.wavefront import ZernikeOPD
    o = detection_optic(d_water, det, t_glass)
    z = ZernikeOPD(o, field=(0.0, 0.0), wavelength=LAM_EM, num_rings=24, zernike_type="fringe", num_terms=37)
    c = np.asarray(z.coeffs, float)
    return {j: float(c[j - 1]) for j in FRINGE_RADIAL}


def best_focus_rms(coeffs: Dict[int, float]) -> float:
    """RMS wavefront (waves) left after refocusing: the symmetric terms above defocus."""
    return math.sqrt(sum(coeffs[j] ** 2 / (FRINGE_RADIAL[j][0] + 1) for j in FRINGE_RADIAL if j != 4))


def w040_analytic_waves(d_water: float, na: float, t_glass: float = 0.17) -> float:
    """Primary spherical aberration of a flat water layer plus window, in waves."""
    w = d_water / 8.0 * (N_WATER_EM ** 2 - 1) / N_WATER_EM ** 3 * na ** 4
    w += t_glass / 8.0 * (N_COVER ** 2 - 1) / N_COVER ** 3 * na ** 4
    return w / (LAM_EM * 1e-3)


def _pupil_phase(coeffs: Dict[int, float], rho: np.ndarray) -> np.ndarray:
    return sum(coeffs[j] * np.polyval(FRINGE_RADIAL[j][1], rho) for j in coeffs)


def psf_profiles(coeffs: Dict[int, float], na: float, n_rho: int = 600) -> Dict[str, object]:
    """Scalar detection PSF in water from the pupil phase: axial scan to the peak, then the
    lateral profile there. Intensities relative to the perfect system's peak."""
    from scipy.special import j0
    rho = (np.arange(n_rho) + 0.5) / n_rho
    drho = 1.0 / n_rho
    phase = 2 * np.pi * _pupil_phase(coeffs, rho)
    cos_t = np.sqrt(1.0 - (na * rho / N_WATER_EM) ** 2)

    def on_axis(delta_um):
        dz = 2 * np.pi * N_WATER_EM * delta_um * (cos_t - 1.0) / LAM_EM
        return abs(2 * np.sum(np.exp(1j * (phase + dz)) * rho) * drho) ** 2
    dz_ideal = 2 * LAM_EM / na ** 2 * N_WATER_EM
    lsa = 24.0 * abs(coeffs.get(9, 0.0)) * LAM_EM * N_WATER_EM / na ** 2      # longitudinal spread of the focus
    half = 6 * dz_ideal + 1.5 * lsa
    deltas = np.linspace(-half, half, int(2 * half / (dz_ideal / 25)) | 1)
    axial = np.array([on_axis(d) for d in deltas])
    k = int(np.argmax(axial))
    best = float(deltas[k])
    r = np.linspace(0, 15 * LAM_EM / na, 1500)
    dzb = 2 * np.pi * N_WATER_EM * best * (cos_t - 1.0) / LAM_EM
    pupil = np.exp(1j * (phase + dzb)) * rho * drho
    lateral = np.array([abs(2 * np.sum(pupil * j0(2 * np.pi * na * ri * rho / LAM_EM))) ** 2 for ri in r])
    norm = np.trapezoid(lateral * r, r)
    mtf = {p: float(np.trapezoid(lateral * j0(2 * np.pi * r / p) * r, r) / norm) for p in (2.0, 4.0)}
    return {"deltas_um": deltas - best, "axial": axial, "r_um": r, "lateral": lateral, "strehl": float(axial[k]), "mtf": mtf}


def fwhm(x: np.ndarray, y: np.ndarray, centred: bool = True) -> float:
    """Full width at half maximum of a peak (centred=False: profile from the axis, r >= 0)."""
    y = np.asarray(y, float)
    k = int(np.argmax(y))
    half = y[k] / 2.0
    idx = np.where(y >= half)[0]
    lo, hi = idx[0], idx[-1]

    def cross(i0, i1):
        return x[i0] + (half - y[i0]) * (x[i1] - x[i0]) / (y[i1] - y[i0])
    right = cross(hi, hi + 1) if hi + 1 < len(y) else x[hi]
    if not centred:
        return 2.0 * right
    left = cross(lo, lo - 1) if lo > 0 else x[lo]
    return float(right - left)


@dataclass
class DetectionQuality:
    d_water: float
    na: float
    rms_waves: float
    strehl: float
    w040_optiland_waves: float
    w040_analytic_waves: float
    lateral_fwhm_um: float
    axial_fwhm_um: float
    mtf_2um: float = 0.0          # modulation of a 2 um period pattern (lateral)
    mtf_4um: float = 0.0
    profiles: Dict[str, object] = field(repr=False, default=None)


def detection_quality(d_water: float, det: Detection, t_glass: float = 0.17) -> DetectionQuality:
    c = zernike_symmetric(d_water, det, t_glass)
    prof = psf_profiles(c, det.na)
    return DetectionQuality(d_water, det.na, best_focus_rms(c), prof["strehl"], 6.0 * abs(c[9]),
                            w040_analytic_waves(d_water, det.na, t_glass),
                            fwhm(prof["r_um"], prof["lateral"], centred=False), fwhm(prof["deltas_um"], prof["axial"]),
                            prof["mtf"][2.0], prof["mtf"][4.0], prof)


def system_axial_fwhm(dq: DetectionQuality, sheet_w_um: float) -> float:
    """Axial FWHM of the detection PSF times the sheet profile (1/e^2 radius sheet_w_um)."""
    d = np.asarray(dq.profiles["deltas_um"])
    prod = np.asarray(dq.profiles["axial"]) * np.exp(-2.0 * d ** 2 / sheet_w_um ** 2)
    return fwhm(d, prod)


def objective_nose_mm(d_water: float, det: Detection, t_glass: float = 0.17) -> Tuple[float, float]:
    """(air gap nose -> window, nose -> sheet) for a focus d_water behind the window (paraxial)."""
    gap = det.wd_air - t_glass / N_COVER - d_water / N_WATER_EM
    return gap, d_water + t_glass + gap


# ----------------------------------------------------------------------------- validation
def validate(verbose: bool = True) -> bool:
    ok = True

    def check(name, got, want, tol, unit=""):
        nonlocal ok
        good = abs(got - want) <= tol
        ok &= good
        if verbose:
            print(f"  {'OK  ' if good else 'FAIL'} {name}: {got:.4g} (expected {want:.4g} +- {tol:g}) {unit}")

    print("illumination: Optiland real rays vs the Gaussian (ABCD) model")
    for key, v in variants().items():
        sh = sheet(v)
        tr = trace_illumination(v)
        check(f"{key} waist on the sheet plane (ABCD)", sh.waist_offset_mm, 0.0, 1e-3, "mm")
        check(f"{key} real-ray focus of the thickness fan vs the Gaussian waist", tr["focus_mm"], sh.waist_offset_mm, 0.3, "mm")
        check(f"{key} sheet half-height, real rays vs Gaussian", tr["half_height_mm"], sh.height_mm / 2, 0.02 * sh.height_mm / 2, "mm")
        check(f"{key} geometric rms spot at focus / Gaussian waist (aberration of the lenses)", tr["rms_spot_um"] / sh.w0_um, 0.0, 0.5)
    print("detection: Optiland wavefront through window + water vs the analytic primary spherical aberration")
    for d, na in ((9.08, 0.10), (3.0, 0.15), (9.08, 0.25)):
        det = Detection(na=na)
        c = zernike_symmetric(d, det)
        tol = 0.05 if na <= 0.15 else 0.12
        check(f"W040 at d {d} mm, NA {na}: Optiland 6 x Z9 / analytic", 6.0 * abs(c[9]) / w040_analytic_waves(d, na), 1.0, tol)
    c0 = zernike_symmetric(0.0, Detection(), t_glass=0.0)
    check("perfect objective, no window, no water: residual rms", best_focus_rms(c0), 0.0, 0.002, "waves")
    o = detection_optic(9.08, Detection())
    o.trace_generic(0.0, 0.0, np.array([0.0]), np.array([1.0]), LAM_EM)
    M, N = float(np.asarray(o.surfaces.M)[2][0]), float(np.asarray(o.surfaces.N)[2][0])
    check("detection NA in air after the stop iteration", abs(M) / math.hypot(M, N), 0.25, 0.002)
    q0 = detection_quality(0.0, Detection(), t_glass=0.0)
    check("perfect NA 0.25: lateral FWHM vs 0.514 lambda/NA", q0.lateral_fwhm_um, 0.514 * LAM_EM / 0.25, 0.03, "um")
    ax = 0.88 * LAM_EM / (N_WATER_EM - math.sqrt(N_WATER_EM ** 2 - 0.25 ** 2))
    check("perfect NA 0.25: axial FWHM vs 0.88 lambda/(n - sqrt(n^2 - NA^2))", q0.axial_fwhm_um, ax, 0.08 * ax, "um")
    print("ALL OK" if ok else "SOME CHECKS FAILED")
    return ok


# ----------------------------------------------------------------------------- numbers + figures
FOV_HALF = None   # filled from Detection()


def report_tables(det: Detection = Detection()) -> Tuple[str, Dict]:
    fx, fy = det.fov_mm
    xe_short, xe_long = fy / 2.0, fx / 2.0
    vs = variants()
    vs45 = variants(45.0)
    vs["K45"] = replace(vs45["K"], key="K45", title=vs45["K"].title)
    data: Dict[str, object] = {"fov_mm": [fx, fy], "variants": {}, "detection": [], "camera": [], "system": []}
    lines = ["## Light sheet per arrangement", "",
             f"Camera field {fx:.2f} x {fy:.2f} mm (10x f{det.f_obj:g} + tube lens f{det.f_tl:g} = {det.mag:.1f}x, "
             f"{det.px_obj_um:.2f} um pixels). Short side along the sheet (x edge {xe_short:.2f} mm), long side along "
             f"the sheet height (z edge {xe_long:.2f} mm).", "",
             "| arrangement | beam at the lens | sheet NA (air) | waist FWHM | length 2 zR (water) | FWHM at x = 0.49 mm | FWHM at x = 0.74 mm | height (1/e2) | brightness at z = 0.74 mm | lens moved from its cell centre |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for key in ("K", "K45", "A", "A0", "A2", "B"):
        v = vs[key]
        sh = sheet(v)
        shift = ", ".join(f"{n.split(',')[0]} {d:+.1f} mm" for n, d in sh.lens_shift_mm.items()) or "-"
        lines.append(f"| {key}: {v.title} | {sh.beam_d_mm:.2f} mm | {sh.na_air:.3f} | {sh.fwhm0_um:.1f} um | {2 * sh.zr_mm:.2f} mm | "
                     f"{float(sh.fwhm_at(xe_short)):.1f} um | {float(sh.fwhm_at(xe_long)):.1f} um | {sh.height_mm:.2f} mm | "
                     f"{100 * sh.height_fraction(xe_long):.0f} % | {shift} |")
        data["variants"][key] = {"title": v.title, "beam_d_mm": sh.beam_d_mm, "na_air": sh.na_air, "w0_um": sh.w0_um,
                                 "fwhm0_um": sh.fwhm0_um, "zr_mm": sh.zr_mm, "height_mm": sh.height_mm,
                                 "fwhm_short_edge_um": float(sh.fwhm_at(xe_short)), "fwhm_long_edge_um": float(sh.fwhm_at(xe_long)),
                                 "height_fraction_long_edge": sh.height_fraction(xe_long),
                                 "irradiance_w_cm2_per_mw": sh.peak_irradiance_w_cm2_per_mw(),
                                 "stations": [{"s": st.s, "lens": st.lens.name, "part": st.lens.part, "cell": st.cell} for st in sh.stations],
                                 "lens_shift_mm": sh.lens_shift_mm, "fibre_cell": v.fibre_cell, "collimator": v.collimator.name,
                                 "s_collimator": v.s_collimator}
    lines += ["", "## Detection through the aquarium window (air objective, perfect lens otherwise)", "",
              "Modulation = contrast with which a lateral grating of that period is imaged (perfect lens = upper bound).", "",
              "| water path window -> sheet | NA | rms wavefront after refocus | Strehl | modulation, 2 um period | modulation, 4 um period | axial FWHM | objective nose -> sheet |",
              "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    dq_cache = {}
    for d in (9.08, 4.08, 3.0, 2.0, 0.0):
        for na in (0.20, 0.25, 0.30):
            tg = 0.17 if d > 0 else 0.0
            q = detection_quality(d, Detection(na=na), tg)
            dq_cache[(d, na)] = q
            gap, nose = objective_nose_mm(d, det, tg)
            label = {9.08: "9.1 mm (today)", 4.08: "4.1 mm", 3.0: "3.0 mm", 2.0: "2.0 mm", 0.0: "0 (water-dipping objective)"}[d]
            lines.append(f"| {label} | {na:.2f} | {q.rms_waves:.3f} waves | {q.strehl:.2f} | {q.mtf_2um:.2f} | {q.mtf_4um:.2f} | "
                         f"{q.axial_fwhm_um:.1f} um | " + (f"{nose:.1f} mm |" if d > 0 else "- |"))
            data["detection"].append({"d_water": d, "na": na, "rms_waves": q.rms_waves, "strehl": q.strehl,
                                      "mtf_2um": q.mtf_2um, "mtf_4um": q.mtf_4um,
                                      "lateral_fwhm_um": q.lateral_fwhm_um, "axial_fwhm_um": q.axial_fwhm_um,
                                      "w040_optiland": q.w040_optiland_waves, "w040_analytic": q.w040_analytic_waves})
    lines += ["", "## Axial resolution of the whole system (detection PSF x sheet), NA 0.25", "",
              "| arrangement | water path | at the sheet waist | at x = 0.49 mm | at x = 0.74 mm |", "| --- | --- | --- | --- | --- |"]
    for key in ("K", "A", "A2", "B"):
        sh = sheet(vs[key])
        for d in (9.08, 3.0):
            q = dq_cache[(d, 0.25)]
            vals = [system_axial_fwhm(q, float(sh.w_at(x))) for x in (0.0, xe_short, xe_long)]
            lines.append(f"| {key} | {d:.1f} mm | {vals[0]:.1f} um | {vals[1]:.1f} um | {vals[2]:.1f} um |")
            data["system"].append({"variant": key, "d_water": d, "axial_fwhm_um": vals})
    lines += ["", "## Camera sampling (10x, NA 0.25, 520 nm: Nyquist pixel = lambda / 4 NA = "
              f"{LAM_EM / (4 * 0.25):.2f} um)", "",
              "| sensor | pixels | pixel | tube lens | pixel in the sample | field | sampling |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for name, px, nx, ny in (("IMX178 (HIK MV-CE060-10UM-PRO, today)", 2.4, 3072, 2048), ("IMX264 / IMX250", 3.45, 2448, 2048),
                             ("IMX174 / IMX249", 5.86, 1920, 1200), ("IMX585", 2.9, 3840, 2160)):
        for ftl in (100.0, 200.0):
            dd = replace(det, pixel_um=px, nx=nx, ny=ny, f_tl=ftl)
            nyq = LAM_EM / (4 * 0.25)
            ratio = dd.px_obj_um / nyq
            samp = "Nyquist or finer" if ratio <= 1.0 else f"{ratio:.1f}x coarser than Nyquist"
            lines.append(f"| {name} | {nx} x {ny} | {px} um | f{ftl:g} | {dd.px_obj_um:.2f} um | {dd.fov_mm[0]:.2f} x {dd.fov_mm[1]:.2f} mm | {samp} |")
            data["camera"].append({"sensor": name, "pixel_um": px, "f_tl": ftl, "px_obj_um": dd.px_obj_um, "fov_mm": dd.fov_mm})
    # infinity space
    h_chief = math.hypot(*det.fov_mm) / 2.0 / det.f_obj * det.obj_to_tl
    r_bundle = det.f_obj * det.na
    lines += ["", f"Infinity space: objective -> tube lens {det.obj_to_tl:g} mm; corner chief ray {h_chief:.1f} mm off axis + "
              f"bundle radius {r_bundle:.1f} mm = {h_chief + r_bundle:.1f} mm vs {det.tl_clear_r:g} mm clear radius "
              f"({'no vignetting' if h_chief + r_bundle <= det.tl_clear_r else 'VIGNETTED at the corners'})."]
    data["infinity_space"] = {"chief_mm": h_chief, "bundle_r_mm": r_bundle}
    return "\n".join(lines) + "\n", data


SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # categorical slots 1-4 (validated palette)
SURFACE, INK, INK2, GRID_C = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"


def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID_C, lw=1)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8, length=0)


def figures(out: Path, det: Detection = Detection()) -> List[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    fx, fy = det.fov_mm
    vs = variants()
    # 1. sheet thickness along the propagation
    fig, ax = plt.subplots(figsize=(10, 4.8), facecolor=SURFACE)
    _style(ax)
    x = np.linspace(-0.8, 0.8, 641)
    ax_ideal = 0.88 * LAM_EM / (N_WATER_EM - math.sqrt(N_WATER_EM ** 2 - det.na ** 2))
    ax.axhline(ax_ideal, color=INK2, lw=1, ls="-", alpha=0.7)
    ax.text(-0.79, ax_ideal + 0.6, f"detection depth of field alone ({ax_ideal:.0f} um, NA {det.na:g})", fontsize=8, color=INK2)
    for xe, lab in ((fy / 2, "short side of the field"), (fx / 2, "long side")):
        for sgn in (-1, 1):
            ax.axvline(sgn * xe, color=INK2, lw=1, alpha=0.35)
        ax.text(xe + 0.01, 36.5, lab, fontsize=7.5, color=INK2, rotation=90, va="top")
    labels = {"K": "K - kit: 4x objective (f 31)", "A": "A and B - same thickness (B is taller)", "A2": "A2 - f11.2 collimator + f50 cylinder"}
    for col, key in zip(SERIES, ("K", "A", "A2")):
        sh = sheet(vs[key])
        y = sh.fwhm_at(x)
        ax.plot(x, y, color=col, lw=2, label=labels[key])
        tag = {"K": "K", "A": "A, B", "A2": "A2"}[key]
        yl = float(y[-1]) if y[-1] < 37 else 36.0
        ax.annotate(tag, (0.8, yl), xytext=(4, 2 if key == "A" else -6), textcoords="offset points", fontsize=8.5, color=INK, va="center")
    ax.set_ylim(0, 38)
    ax.set_xlim(-0.8, 0.8)
    ax.set_xlabel("position along the sheet, from the waist (mm)", fontsize=9, color=INK2)
    ax.set_ylabel("sheet thickness FWHM (um)", fontsize=9, color=INK2)
    leg = ax.legend(loc="upper center", frameon=False, fontsize=8, ncol=2)
    for t in leg.get_texts():
        t.set_color(INK)
    fig.suptitle("The kit's sheet is thin but only 0.07 mm long; A, A2 and B stay below 13 um over the field's short side",
                 fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 0.97, 0.95))
    p = out / "ls_sheet_thickness.png"
    fig.savefig(p, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    paths.append(p)
    # 2. detection Strehl vs water path
    fig, ax = plt.subplots(figsize=(10, 4.6), facecolor=SURFACE)
    _style(ax)
    ds = np.linspace(0.25, 12.0, 48)
    for col, na in zip(SERIES, (0.20, 0.25, 0.30)):
        st = [detection_quality(d, Detection(na=na)).strehl for d in ds]
        ax.plot(ds, st, color=col, lw=2, label=f"NA {na:.2f}")
        ax.annotate(f"NA {na:.2f}", (ds[-1], st[-1]), xytext=(4, 0), textcoords="offset points", fontsize=8.5, color=INK, va="center")
    for d, lab in ((9.08, "today: 9.1 mm"), (3.0, "proposed: 3 mm")):
        ax.axvline(d, color=INK2, lw=1, alpha=0.5)
        ax.text(d + 0.1, 0.03, lab, fontsize=8, color=INK2, va="bottom")
    ax.axhline(0.8, color=INK2, lw=1, ls=":", alpha=0.7)
    ax.text(0.3, 0.81, "diffraction limited (Strehl 0.8)", fontsize=8, color=INK2)
    ax.set_xlim(0, 13.2)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("water between the detection window and the sheet (mm)", fontsize=9, color=INK2)
    ax.set_ylabel("Strehl ratio at best focus", fontsize=9, color=INK2)
    leg = ax.legend(loc="upper right", frameon=False, fontsize=8)
    for t in leg.get_texts():
        t.set_color(INK)
    fig.suptitle("An air objective looking through 9 mm of water loses most of its sharpness; 3 mm keeps NA 0.25 near the limit",
                 fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 0.97, 0.95))
    p = out / "ls_detection_strehl.png"
    fig.savefig(p, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    paths.append(p)
    # 3. detection PSFs
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE)
    for ax in axs:
        _style(ax)
    for col, (d, tg, lab) in zip(SERIES, ((0.0, 0.0, "no water (dipping objective)"), (3.0, 0.17, "3 mm water"), (9.08, 0.17, "9.1 mm water (today)"))):
        q = detection_quality(d, det, tg)
        axs[0].plot(np.r_[-q.profiles["r_um"][::-1], q.profiles["r_um"]], np.r_[q.profiles["lateral"][::-1], q.profiles["lateral"]], color=col, lw=2, label=lab)
        axs[1].plot(q.profiles["deltas_um"], q.profiles["axial"], color=col, lw=2, label=lab)
    axs[0].set_xlim(-4, 4)
    axs[1].set_xlim(-80, 80)
    axs[0].set_xlabel("lateral position (um)", fontsize=9, color=INK2)
    axs[1].set_xlabel("axial position (um)", fontsize=9, color=INK2)
    axs[0].set_ylabel("intensity (perfect system = 1)", fontsize=9, color=INK2)
    leg = axs[1].legend(loc="upper right", frameon=False, fontsize=8)
    for t in leg.get_texts():
        t.set_color(INK)
    fig.suptitle(f"Detection PSF of the 10x (NA {det.na:g}) through the window and water, at best focus", fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = out / "ls_detection_psf.png"
    fig.savefig(p, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    paths.append(p)
    return paths


# ----------------------------------------------------------------------------- export
def _fmt(v) -> str:
    if isinstance(v, str):
        if v.startswith("ideal:"):
            return f"IdealMaterial(n={float(v.split(':')[1])})"
        return '"' + v + '"'
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):
        return "np.inf" if math.isinf(v) else f"{v:.6f}"
    return str(v)


def script(specs, aperture, lam, title: str, fields=(0.0,), field_type="angle", fan: str = "y") -> str:
    body = "\n".join(f"lens.surfaces.add(index={i}, " + ", ".join(f"{k}={_fmt(v)}" for k, v in kw.items()) + ")" for i, kw in specs)
    flds = "\n".join(f"lens.fields.add(y={f})" for f in fields)
    q = chr(34) * 3
    return (f"{q}\n{title}\n\nExported by openuc2-lightsheet/ls_optics.py. Unfolded path (the fold mirror is left out), "
            f"lengths in mm.\n{q}\nimport numpy as np\nfrom optiland import optic\nfrom optiland.materials import IdealMaterial\n\n"
            f"lens = optic.Optic()\n{body}\nlens.set_aperture(aperture_type=\"{aperture[0]}\", value={aperture[1]:.6f})\n"
            f"lens.fields.set_type(\"{field_type}\")\n{flds}\nlens.wavelengths.add(value={lam}, is_primary=True)\n\n"
            "if __name__ == \"__main__\":\n    import matplotlib\n    matplotlib.use(\"Agg\")\n    import matplotlib.pyplot as plt\n"
            f"    res = lens.draw(num_rays=9, projection=\"{'YZ' if fan == 'y' else 'XZ'}\")\n"
            "    fig = res[0] if isinstance(res, tuple) else plt.gcf()\n"
            "    fig.savefig(__file__.replace('.py', '.png'), dpi=150, bbox_inches='tight')\n"
            "    print('saved', __file__.replace('.py', '.png'))\n")


def export(out: Path) -> List[Path]:
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for key, v in variants().items():
        d = out / key
        d.mkdir(parents=True, exist_ok=True)
        specs = illumination_specs(v)
        for fan in ("y", "z"):
            p = d / f"ls_illumination_{key}_{'thickness' if fan == 'y' else 'height'}_optiland.py"
            p.write_text(script(specs, ("objectNA", math.sin(v.fibre.theta())), LAM_EX,
                                f"openUC2 light sheet, arrangement {key}: {v.title}\nIllumination, fibre to the sheet plane; "
                                f"drawn in the {'sheet-thickness (YZ)' if fan == 'y' else 'sheet-height (XZ)'} plane.", fan=fan), encoding="utf-8")
            paths.append(p)
    for d in (9.08, 3.0):
        p = out / f"ls_detection_water{d:g}mm_optiland.py"
        o = detection_optic(d, Detection())
        p.write_text(script(detection_specs(d, Detection()), ("float_by_stop_size", float(o.aperture.value)), LAM_EM,
                            f"openUC2 light sheet: detection through {d:g} mm of water and a 0.17 mm window with a perfect 10x "
                            f"(f 20, NA 0.25) and the f100 tube lens.", fields=(0.0, 0.7), field_type="object_height"), encoding="utf-8")
        paths.append(p)
    return paths


def main() -> None:
    ap = argparse.ArgumentParser(description="openUC2 light-sheet optics")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--report", default=None, metavar="FILE.md")
    ap.add_argument("--json", default=None, metavar="FILE.json")
    ap.add_argument("--figures", default=None, metavar="DIR")
    ap.add_argument("--export", default=None, metavar="DIR")
    a = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if a.validate:
        if not validate():
            raise SystemExit(1)
    if a.report or a.json:
        text, data = report_tables()
        if a.report:
            Path(a.report).parent.mkdir(parents=True, exist_ok=True)
            Path(a.report).write_text(text, encoding="utf-8")
        if a.json:
            Path(a.json).write_text(json.dumps(data, indent=2, default=lambda o: list(o) if isinstance(o, tuple) else str(o)), encoding="utf-8")
        print(text)
    if a.figures:
        for p in figures(Path(a.figures)):
            print("wrote", p)
    if a.export:
        for p in export(Path(a.export)):
            print("wrote", p)


if __name__ == "__main__":
    main()
