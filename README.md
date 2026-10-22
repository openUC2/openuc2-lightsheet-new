# openuc2-lightsheet

Optics and arrangement of the openUC2 light-sheet microscope (LightSheet add-on with the FluoBox
laser, CoreBox, InfinityBox), ahead of putting it into cubes. Round 1 = simulation and module
arrangement; design note: [`docs/lightsheet.md`](docs/lightsheet.md).

| file | what it does |
| --- | --- |
| `ls_optics.py` | the sheet as a Gaussian beam per axis, Optiland models of illumination and detection (water + window), validation, tables, figures, stand-alone Optiland scripts |
| `ls_layout.py` | the arrangements K (today), A, A2, B on the 50 mm grid: drawings in the style of the kit sketch, cell-by-cell module lists with CAD-new names |
| `docs/lightsheet.md` | findings, what to change, samples and mounts, open questions |
| `docs/ls_numbers.md` / `.json` | all numbers |
| `layouts/<key>/` | `ls_layout.png/json`, exported Optiland scripts |

```bash
uv run --with optiland --with numpy --with scipy --with matplotlib python ls_optics.py --validate
uv run --with optiland --with numpy --with scipy --with matplotlib python ls_optics.py --report docs/ls_numbers.md --json docs/ls_numbers.json --figures docs --export layouts
uv run --with optiland --with numpy --with scipy --with matplotlib python ls_layout.py
```
