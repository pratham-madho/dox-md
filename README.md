# DOX⁺–membrane molecular dynamics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24459737.svg)](https://doi.org/10.5281/zenodo.24459737)

GROMACS molecular dynamics simulations of protonated doxorubicin (DOX⁺) interacting with lipid bilayers:

- **System A**: fluid DPPC bilayer (`system_A/`)
- **System B**: DPPC/cholesterol bilayer (`system_B/`)

The write-up is in [`report/report.pdf`](report/report.pdf) (source: `report/report.tex`, figures from `report/make_figures.py`).

## Layout

| Path | Contents |
|---|---|
| `system_A/`, `system_B/` | GROMACS inputs, run outputs, analysis scripts and `.xvg` data |
| `charmm-gui-*/` , `charmm-gui.tgz` | CHARMM-GUI Membrane Builder / Ligand Reader outputs used to build the systems |
| `DOX_plus.sdf`, `Conformer3D_COMPOUND_CID_31703.sdf`, `protonate.py` | Doxorubicin structure (PubChem CID 31703) and protonation script |
| `report/` | Report, figures and figure script |

## Trajectories

Production trajectories (`*.xtc`, ~7 GB) and rendered videos (`*.mp4`) are too large for GitHub and are archived separately:

> **Data:** [doi:10.5281/zenodo.24459737](https://doi.org/10.5281/zenodo.24459737) (Zenodo)

## Environment

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Python 3.12; main packages: MDAnalysis, RDKit, NumPy, SciPy, Matplotlib.

## References

- Reference paper: [arXiv:2109.00436](https://arxiv.org/abs/2109.00436)
- CHARMM-GUI: https://www.charmm-gui.org
