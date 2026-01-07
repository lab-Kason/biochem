"""
Simple Molecular Docking Example Using RDKit + OpenMM
No need to install AutoDock Vina - uses packages you already have!

This demonstrates conformer generation and energy-based ranking,
which is conceptually similar to what Vina does.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np
import pandas as pd

print("🧬 Simple Docking with RDKit (No Vina Needed!)")
print("=" * 60)

# 1. Create a ligand from SMILES (or load from file)
print("\n1️⃣ Creating ligand molecule...")
smiles = 'CC(=O)Oc1ccccc1C(=O)O'  # Aspirin
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

print(f"   Molecule: Aspirin")
print(f"   SMILES: {smiles}")
print(f"   Atoms: {mol.GetNumAtoms()}")
print(f"   Molecular Weight: {Descriptors.MolWt(mol):.2f}")

# 2. Generate multiple conformers (like Vina generates poses)
print("\n2️⃣ Generating conformers (like docking poses)...")
num_conformers = 10
AllChem.EmbedMultipleConfs(
    mol, 
    numConfs=num_conformers,
    randomSeed=42,
    useRandomCoords=True
)
print(f"   Generated {mol.GetNumConformers()} conformers")

# 3. Optimize each conformer with MMFF force field
print("\n3️⃣ Optimizing conformers...")
energies = []
for conf_id in range(mol.GetNumConformers()):
    # Minimize energy
    result = AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
    
    # Calculate energy
    props = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    energy = ff.CalcEnergy()
    energies.append(energy)
    
    print(f"   Conformer {conf_id + 1}: {energy:.2f} kcal/mol")

# 4. Rank by energy (like Vina ranks by binding affinity)
print("\n4️⃣ Ranking conformers by energy...")
energies = np.array(energies)
sorted_indices = np.argsort(energies)

results = pd.DataFrame({
    'Conformer': sorted_indices + 1,
    'Energy (kcal/mol)': energies[sorted_indices],
    'Rank': range(1, len(energies) + 1)
})

print(results.to_string(index=False))

# 5. Save best conformer
best_conf_id = sorted_indices[0]
best_energy = energies[best_conf_id]

print(f"\n5️⃣ Best conformer: #{best_conf_id + 1}")
print(f"   Energy: {best_energy:.2f} kcal/mol")

# Save as PDB
Chem.MolToPDBFile(mol, 'best_conformer.pdb', confId=int(best_conf_id))
print(f"   Saved to: best_conformer.pdb")

# Save all conformers as SDF
writer = Chem.SDWriter('all_conformers.sdf')
for conf_id in range(mol.GetNumConformers()):
    writer.write(mol, confId=conf_id)
writer.close()
print(f"   All conformers saved to: all_conformers.sdf")

# 6. Visualization hint
print("\n6️⃣ Visualize results:")
print("   • In notebook: import py3Dmol")
print("   • Command line: pymol best_conformer.pdb")
print("   • Online: https://molstar.org/viewer/")

print("\n" + "=" * 60)
print("✅ Docking simulation complete!")
print("=" * 60)
print("\n📝 Notes:")
print("   • This uses force field energies (not binding affinity)")
print("   • For actual docking, see VINA_ALTERNATIVES.md")
print("   • GNINA recommended for real projects")
print("\n💡 But this approach works for:")
print("   • Understanding conformational flexibility")
print("   • Generating starting structures")
print("   • Quick energy estimates")
print("   • Educational purposes")
