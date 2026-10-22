"""
openUC2 light sheet, arrangement A: compact 3 x 4: f28 collimator, 2:1 Galilean reducer (f100, f-50), f50 cylinder
Illumination, fibre to the sheet plane; drawn in the sheet-thickness (YZ) plane.

Exported by openuc2-lightsheet/ls_optics.py. Unfolded path (the fold mirror is left out), lengths in mm.
"""
import numpy as np
from optiland import optic
from optiland.materials import IdealMaterial

lens = optic.Optic()
lens.surfaces.add(index=0, radius=np.inf, thickness=28.000000, comment="fibre tip")
lens.surfaces.add(index=1, surface_type="paraxial", f=28.000000, comment="collimator f28 D16 (FluoBox laser)", thickness=50.000000, is_stop=True)
lens.surfaces.add(index=2, surface_type="paraxial", f=100.000000, comment="lens f100 D40 (CoreBox)", thickness=50.000000)
lens.surfaces.add(index=3, surface_type="paraxial", f=-50.000000, comment="lens f-50 D43 (CoreBox)", thickness=95.430302)
lens.surfaces.add(index=4, surface_type="toroidal", radius_x=np.inf, radius_y=26.120000, material="N-BK7", comment="cylinder f50, power across the sheet, curved face", thickness=5.000000)
lens.surfaces.add(index=5, radius=np.inf, comment="cylinder f50, power across the sheet, flat face", thickness=38.319698)
lens.surfaces.add(index=6, radius=np.inf, material=IdealMaterial(n=1.523), comment="window, outer face", thickness=0.170000)
lens.surfaces.add(index=7, radius=np.inf, material=IdealMaterial(n=1.3371), comment="window, inner face (water)", thickness=11.080000)
lens.surfaces.add(index=8, radius=np.inf, comment="sheet plane (aquarium cube centre)")
lens.set_aperture(aperture_type="objectNA", value=0.086190)
lens.fields.set_type("angle")
lens.fields.add(y=0.0)
lens.wavelengths.add(value=0.488, is_primary=True)

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res = lens.draw(num_rays=9, projection="YZ")
    fig = res[0] if isinstance(res, tuple) else plt.gcf()
    fig.savefig(__file__.replace('.py', '.png'), dpi=150, bbox_inches='tight')
    print('saved', __file__.replace('.py', '.png'))
