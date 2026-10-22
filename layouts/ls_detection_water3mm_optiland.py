"""
openUC2 light sheet: detection through 3 mm of water and a 0.17 mm window with a perfect 10x (f 20, NA 0.25) and the f100 tube lens.

Exported by openuc2-lightsheet/ls_optics.py. Unfolded path (the fold mirror is left out), lengths in mm.
"""
import numpy as np
from optiland import optic
from optiland.materials import IdealMaterial

lens = optic.Optic()
lens.surfaces.add(index=0, radius=np.inf, thickness=3.000000, material=IdealMaterial(n=1.3354), comment="sample point in water")
lens.surfaces.add(index=1, radius=np.inf, thickness=0.170000, material=IdealMaterial(n=1.523), comment="detection window")
lens.surfaces.add(index=2, radius=np.inf, thickness=17.641860, comment="air gap to the objective")
lens.surfaces.add(index=3, surface_type="paraxial", f=20.000000, thickness=60.000000, is_stop=True, comment="10x objective (perfect)")
lens.surfaces.add(index=4, surface_type="paraxial", f=100.000000, thickness=100.000000, comment="tube lens f100")
lens.surfaces.add(index=5, radius=np.inf, comment="camera")
lens.set_aperture(aperture_type="float_by_stop_size", value=10.180166)
lens.fields.set_type("object_height")
lens.fields.add(y=0.0)
lens.fields.add(y=0.7)
lens.wavelengths.add(value=0.52, is_primary=True)

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res = lens.draw(num_rays=9, projection="YZ")
    fig = res[0] if isinstance(res, tuple) else plt.gcf()
    fig.savefig(__file__.replace('.py', '.png'), dpi=150, bbox_inches='tight')
    print('saved', __file__.replace('.py', '.png'))
