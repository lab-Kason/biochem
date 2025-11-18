Getting started: Python for docking and molecular dynamics

Overview
- This project provides two small example scripts showing how Python can be used to run a docking workflow (wrapper around AutoDock Vina) and a minimal molecular dynamics example (OpenMM).
- Two installation paths are supported: recommended `conda` (best for RDKit/OpenMM) and `pip` (lighter; may not work for some packages on macOS). If you don't have conda, I recommend installing Miniconda from https://docs.conda.io/en/latest/miniconda.html

Recommended (conda) installation
1. Install Miniconda or Anaconda for macOS.
2. Create an environment and install packages:
```bash
conda create -n md-dock python=3.11 -y
conda activate md-dock
conda install -c conda-forge openmm rdkit mdanalysis mdtraj nglview autodock-vina -y
```

Pip / venv installation (if you cannot use conda)
1. Create and activate a virtualenv:
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```
2. Install Python packages (note: RDKit and OpenMM may be problematic via pip on macOS; prefer conda):
```bash
pip install openmm mdtraj mdanalysis py3Dmol vina
# Optional attempt for RDKit (may fail):
pip install rdkit-pypi || echo "rdkit install failed: use conda instead"
```

Validation
- Check OpenMM: `python -c "import openmm; print(openmm.__version__)"`
- Check Vina: `vina --help` (CLI) or `python -c "import vina; print(vina.__version__)"` (Python wrapper)

Running the examples
- Docking example:
  - Prepare or place a receptor file `receptor.pdbqt` and a ligand (smiles or `ligand.sdf`). The script will try to create a 3D ligand and call Vina if available.
  - Run: `python docking_example.py`
- OpenMM example:
  - Runs a short, self-contained test simulation if `openmm` is installed.
  - Run: `python openmm_example.py`

Notes & troubleshooting
- Many cheminformatics and MD packages are easiest to install with `conda` on macOS. If you hit errors installing with `pip`, install Miniconda and use the recommended conda commands above.
- If `autodock-vina` (CLI) is not installed, consider `conda install -c conda-forge autodock-vina` or download a binary from the AutoDock site.

References
- OpenMM docs: https://openmm.org
- AutoDock Vina: https://vina.scripps.edu
- RDKit: https://www.rdkit.org
