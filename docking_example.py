#!/usr/bin/env python3
"""
Simple docking example (illustrative).

What it does:
- Demonstrates ligand preparation (RDKit) to 3D and attempts to call AutoDock Vina.
- If `vina` Python package is available it uses that; otherwise tries the `vina` CLI.

Notes:
- You need a receptor in `receptor.pdbqt` (prepared by AutoDockTools or OpenBabel).
- The script will not automatically install Vina or prepare the receptor; this is a minimal example.
"""
import os
import shutil
import subprocess
import sys

def prepare_ligand_from_smiles(smiles: str, out_molfile: str = "ligand.sdf"):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as e:
        print("RDKit not available. Install RDKit or prepare ligand file manually.")
        raise

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    w = Chem.SDWriter(out_molfile)
    w.write(mol)
    w.close()
    print(f"Wrote ligand to {out_molfile}")

def try_vina_python(receptor, ligand, center, size, out='out.pdbqt'):
    try:
        from vina import Vina
    except Exception:
        return False
    v = Vina(sf_name='vina')
    v.set_receptor(receptor)
    v.set_ligand_from_file(ligand)
    v.compute_vina_maps(center=center, box_size=size)
    v.dock(n_poses=5)
    v.write_poses(out, n_poses=5)
    print('Wrote poses to', out)
    return True

def try_vina_cli(receptor, ligand, center, size, out='out.pdbqt'):
    exe = None
    if shutil.which('vina'):
        exe = 'vina'
    elif shutil.which('smina'):
        exe = 'smina'
    if exe is None:
        print('Vina CLI not found in PATH (tried `vina` and `smina`).')
        return False
    
    # Convert ligand to PDBQT if needed (vina v1.2.7 requires PDBQT for both receptor and ligand)
    ligand_pdbqt = ligand.replace('.sdf', '.pdbqt')
    if not ligand.endswith('.pdbqt'):
        if shutil.which('obabel'):
            print(f'Converting {ligand} to {ligand_pdbqt}...')
            subprocess.run(['obabel', '-isdf', ligand, '-opdbqt', '-O', ligand_pdbqt, 
                          '--partialcharge', 'gasteiger'], check=True)
            ligand = ligand_pdbqt
        else:
            print('OpenBabel not found - cannot convert ligand to PDBQT')
            return False
    
    cmd = [
        exe,
        '--receptor', receptor,
        '--ligand', ligand,
        '--center_x', str(center[0]),
        '--center_y', str(center[1]),
        '--center_z', str(center[2]),
        '--size_x', str(size[0]),
        '--size_y', str(size[1]),
        '--size_z', str(size[2]),
        '--out', out,
        '--exhaustiveness', '8'
    ]
    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)
    print('Wrote poses to', out)
    return True

def main():
    receptor = 'receptor.pdbqt'
    ligand_smi = 'CCO'  # ethanol, tiny example
    ligand = 'ligand.sdf'
    center = [0,0,0]
    size = [20,20,20]

    # If receptor not present, we'll build a toy receptor later in the script.

    # Prepare ligand
    try:
        prepare_ligand_from_smiles(ligand_smi, ligand)
    except Exception:
        print('Could not prepare ligand automatically. Place a ligand file named', ligand)

    # If no receptor file, build a tiny toy receptor and convert it to PDBQT using OpenBabel
    if not os.path.exists(receptor):
        print(f'{receptor} not found — building a toy receptor `receptor.pdb` and converting to pdbqt (requires openbabel).')
        # build simple receptor: a small cluster of carbon atoms
        with open('receptor.pdb', 'w') as f:
            f.write('REMARK toy receptor\n')
            idx = 1
            for x in [ -1.0, 0.0, 1.0 ]:
                for y in [ -1.0, 0.0, 1.0 ]:
                    for z in [ -1.0, 0.0, 1.0 ]:
                        f.write(f"HETATM{idx:5d}  C   UNK A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
                        idx += 1
        # convert to pdbqt using obabel (openbabel)
        if shutil.which('obabel'):
            cmd = ['obabel', '-ipdb', 'receptor.pdb', '-opdbqt', '-O', 'receptor.pdbqt', '--partialcharge', 'gasteiger']
            print('Running:', ' '.join(cmd))
            subprocess.run(cmd, check=True)
        else:
            print('OpenBabel not found in PATH. Try `conda run -n md-dock obabel ...` or install openbabel.')

    # Try vina python API, then CLI
    if try_vina_python(receptor, ligand, center, size):
        return
    if try_vina_cli(receptor, ligand, center, size):
        return

    print('No Vina backend available. Install AutoDock Vina (CLI) or build/install the `vina` Python package with its C++ dependencies.')

if __name__ == '__main__':
    main()
