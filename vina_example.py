"""
AutoDock Vina Docking Example
Now that Vina is installed, here's how to use it!
"""

import subprocess
import pandas as pd
from pathlib import Path

print("🧬 AutoDock Vina Quick Start Guide")
print("=" * 60)

# Check Vina installation
result = subprocess.run(['vina', '--version'], capture_output=True, text=True)
print(f"\n✅ AutoDock Vina installed: {result.stdout.strip()}")

print("\n" + "=" * 60)
print("📚 Basic Vina Workflow")
print("=" * 60)

print("""
1️⃣ Prepare Receptor (Protein)
   ────────────────────────────
   # Remove water, add hydrogens, save as PDBQT
   prepare_receptor -r protein.pdb -o receptor.pdbqt

2️⃣ Prepare Ligand (Small Molecule)
   ────────────────────────────
   # From SMILES or PDB, convert to PDBQT
   obabel -:"CCO" -O ligand.pdb --gen3d
   prepare_ligand -l ligand.pdb -o ligand.pdbqt

3️⃣ Define Search Space
   ────────────────────────────
   # Find binding site coordinates (center_x, center_y, center_z)
   # Define box size (size_x, size_y, size_z) in Angstroms

4️⃣ Run Docking
   ────────────────────────────
   vina --receptor receptor.pdbqt \\
        --ligand ligand.pdbqt \\
        --center_x 0 --center_y 0 --center_z 0 \\
        --size_x 20 --size_y 20 --size_z 20 \\
        --out docked.pdbqt \\
        --exhaustiveness 8

5️⃣ Analyze Results
   ────────────────────────────
   # Results are in docked.pdbqt
   # View with PyMOL or parse programmatically
""")

print("\n" + "=" * 60)
print("💻 Example: Quick Docking Test")
print("=" * 60)

# Check if we have the example files
receptor_file = Path('receptor.pdbqt')
ligand_file = Path('ligand.pdbqt')

if receptor_file.exists() and ligand_file.exists():
    print("\n✅ Example files found!")
    print(f"   Receptor: {receptor_file}")
    print(f"   Ligand: {ligand_file}")
    
    print("\n🚀 Running docking...")
    print("   Command:")
    cmd = [
        'vina',
        '--receptor', str(receptor_file),
        '--ligand', str(ligand_file),
        '--center_x', '0',
        '--center_y', '0', 
        '--center_z', '0',
        '--size_x', '20',
        '--size_y', '20',
        '--size_z', '20',
        '--out', 'vina_docked.pdbqt',
        '--exhaustiveness', '8'
    ]
    print("   " + " ".join(cmd))
    
    # Run docking
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print("\n✅ Docking complete!")
        print(f"   Output: vina_docked.pdbqt")
        
        # Parse results
        print("\n📊 Docking Results:")
        with open('vina_docked.pdbqt', 'r') as f:
            for line in f:
                if line.startswith('REMARK VINA RESULT:'):
                    parts = line.split()
                    affinity = parts[3]
                    print(f"   Binding affinity: {affinity} kcal/mol")
                    break
    else:
        print("\n❌ Docking failed:")
        print(result.stderr)
else:
    print("\n⚠️  Example files not found.")
    print("   Create receptor.pdbqt and ligand.pdbqt to test docking")

print("\n" + "=" * 60)
print("🔧 Useful Vina Commands")
print("=" * 60)

print("""
# Show help
vina --help

# Quick docking with default settings
vina --receptor protein.pdbqt --ligand ligand.pdbqt \\
     --center_x 0 --center_y 0 --center_z 0 \\
     --size_x 20 --size_y 20 --size_z 20 \\
     --out result.pdbqt

# High-quality docking (slower)
vina --receptor protein.pdbqt --ligand ligand.pdbqt \\
     --center_x 0 --center_y 0 --center_z 0 \\
     --size_x 20 --size_y 20 --size_z 20 \\
     --exhaustiveness 32 \\
     --num_modes 20 \\
     --out result.pdbqt

# Batch docking (multiple ligands)
for lig in ligands/*.pdbqt; do
    vina --receptor protein.pdbqt --ligand "$lig" \\
         --center_x 0 --center_y 0 --center_z 0 \\
         --size_x 20 --size_y 20 --size_z 20 \\
         --out "results/$(basename $lig .pdbqt)_docked.pdbqt"
done
""")

print("\n" + "=" * 60)
print("📖 More Resources")
print("=" * 60)
print("""
• Vina Manual: https://autodock-vina.readthedocs.io/
• Prepare files: https://autodock-vina.readthedocs.io/en/latest/docking_basic.html
• Python API: import vina (if installed via pip)
• Your notebook: visualize_md_docking.ipynb (already set up!)
""")

print("\n✅ You're ready to dock!")
