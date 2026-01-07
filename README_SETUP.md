# 🧬 Complete Molecular Modeling Setup Guide

**Last Updated:** November 19, 2025  
**Environment:** md-dock (Python 3.11, conda)

---

## ✅ What's Installed and Working

### Core Packages (11/11)
- **OpenMM 8.4** - MD simulations (CPU, OpenCL platforms)
- **MDAnalysis 2.10.0** - Trajectory analysis  
- **BioPython 1.86** - PDB files, sequences
- **RDKit 2023.09.6** - Molecule generation, cheminformatics
- **ProLIF 2.0.3** - Protein-ligand interaction analysis
- **NGLView 4.0** - 3D visualization (JupyterLab)
- **Py3Dmol 2.5.3** - 3D visualization (VS Code)
- **NumPy, SciPy, Pandas, Matplotlib** - Scientific computing

### Command-Line Tools
- **AutoDock Vina 1.2.7** - Molecular docking (you installed Nov 17)
- **Open Babel 3.1.0** - File format conversion
- **PDBFixer** - Structure preparation (in base env)

---

## 🚀 Quick Start (3 Commands)

```bash
# 1. Activate environment
conda activate md-dock

# 2. Test everything
bash test_md_dock.sh

# 3. Open your notebook
code visualize_md_docking.ipynb
```

---

## 💻 Common Tasks

### 1. Molecular Docking (AutoDock Vina)
```bash
vina --receptor protein.pdbqt \
     --ligand drug.pdbqt \
     --center_x 0 --center_y 0 --center_z 0 \
     --size_x 20 --size_y 20 --size_z 20 \
     --out result.pdbqt
```

### 2. File Format Conversion (Open Babel)
```bash
# SMILES to 3D
obabel -:"CCO" -O ethanol.pdb --gen3d

# Convert formats
obabel input.sdf -O output.pdbqt
```

### 3. Molecule Generation (RDKit)
```python
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles('CCO')
AllChem.EmbedMolecule(mol)
Chem.MolToPDBFile(mol, 'output.pdb')
```

### 4. MD Simulation (OpenMM)
```python
from openmm.app import *
from openmm import *
from openmm.unit import *

forcefield = ForceField('amber14-all.xml')
system = forcefield.createSystem(topology)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(topology, system, integrator)
simulation.step(10000)
```

### 5. Trajectory Analysis (MDAnalysis)
```python
import MDAnalysis as mda

u = mda.Universe('topology.pdb', 'trajectory.dcd')
# Calculate RMSD, RMSF, distances, etc.
```

### 6. 3D Visualization (Py3Dmol - VS Code friendly)
```python
import py3Dmol

view = py3Dmol.view(width=800, height=600)
view.addModel(pdb_string, 'pdb')
view.setStyle({'cartoon': {}})
view.show()
```

---

## 🔧 Troubleshooting

### "Module not found" errors
**Problem:** Your shell uses Python 3.14, but packages are in Python 3.11  
**Solution:**
```bash
conda activate md-dock  # Always activate first!
```

### NGLView not displaying in VS Code
**Solution:** Use Py3Dmol instead (works better in VS Code)

### Need more packages
```bash
conda activate md-dock
conda install -c conda-forge <package-name>
```

---

## 📚 Key Resources

**Documentation:**
- OpenMM: http://docs.openmm.org/
- Vina: https://autodock-vina.readthedocs.io/
- MDAnalysis: https://userguide.mdanalysis.org/
- RDKit: https://www.rdkit.org/docs/

**Databases:**
- PDB: https://www.rcsb.org/
- ZINC: https://zinc.docking.org/
- PubChem: https://pubchem.ncbi.nlm.nih.gov/

**Your Files:**
- `visualize_md_docking.ipynb` - Main notebook with examples
- `test_md_dock.sh` - Test your installation
- `simple_docking_example.py` - RDKit conformer generation
- `vina_example.py` - AutoDock Vina tutorial

---

## 🎓 Learning Path

**Week 1: Basics**
1. Run `bash test_md_dock.sh`
2. Open `visualize_md_docking.ipynb`
3. Try examples in the notebook

**Week 2: Real Work**
1. Download protein from PDB
2. Prepare for docking
3. Run docking and analyze

**Week 3: Advanced**
1. MD simulations with OpenMM
2. Trajectory analysis
3. Interaction analysis with ProLIF

---

## ⚠️ Important Notes

1. **Always activate conda environment first:** `conda activate md-dock`
2. **Vina installation:** You installed this yourself on Nov 17 (not via conda)
3. **Python version:** Use Python 3.11 in md-dock environment, not system Python 3.14
4. **VS Code kernel:** Select "Python 3.11.0 ('md-dock')" for notebooks

---

## 📊 Installation Status

Run `bash test_md_dock.sh` to see:
```
✅ OpenMM 8.4
✅ MDAnalysis 2.10.0  
✅ BioPython 1.86
✅ RDKit 2023.09.6
✅ ProLIF 2.0.3
✅ NGLView 4.0
✅ Py3Dmol 2.5.3
✅ NumPy, SciPy, Pandas, Matplotlib
✅ AutoDock Vina 1.2.7
✅ Open Babel 3.1.0

📊 Summary: 11/11 packages working
```

---

## ✅ You're Ready!

Your environment is complete. Start with:
```bash
conda activate md-dock
code visualize_md_docking.ipynb
```

**Everything you need for molecular modeling is installed and working!** 🎉
