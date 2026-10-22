"""
openUC2 light sheet, arrangement B: 4 x 4, tall uniform sheet: f28 collimator, f100 cylinder two cells before the sheet
Illumination, fibre to the sheet plane; drawn in the sheet-height (XZ) plane.

Exported by openuc2-lightsheet/ls_optics.py. Unfolded path (the fold mirror is left out), lengths in mm.
"""
import numpy as np
from optiland import optic
from optiland.materials import IdealMaterial

lens = optic.Optic()
lens.surfaces.add(index=0, radius=np.inf, thickness=28.000000, comment="fibre tip")
lens.surfaces.add(index=1, surface_type="paraxial", f=28.000000, comment="collimator f28 D16 (FluoBox laser)", thickness=95.432719, is_stop=True)
lens.surfaces.add(index=2, surface_type="toroidal", radius_x=np.inf, radius_y=52.240000, material="N-BK7", comment="cylinder f100, power across the sheet, curved face", thickness=5.000000)
lens.surfaces.add(index=3, radius=np.inf, comment="cylinder f100, power across the sheet, flat face", thickness=88.317281)
lens.surfaces.add(index=4, radius=np.inf, material=IdealMaterial(n=1.523), comment="window, outer face", thickness=0.170000)
lens.surfaces.add(index=5, radius=np.inf, material=IdealMaterial(n=1.3371), comment="window, inner face (water)", thickness=11.080000)
lens.surfaces.add(index=6, radius=np.inf, comment="sheet plane (aquarium cube centre)")
lens.set_aperture(aperture_type="objectNA", value=0.086190)
lens.fields.set_type("angle")
lens.fields.add(y=0.0)
lens.wavelengths.add(value=0.488, is_primary=True)

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res = lens.draw(num_rays=9, projection="XZ")
    fig = res[0] if isinstance(res, tuple) else plt.gcf()
    fig.savefig(__file__.replace('.py', '.png'), dpi=150, bbox_inches='tight')
    print('saved', __file__.replace('.py', '.png'))
