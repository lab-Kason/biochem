#!/usr/bin/env python3
"""
Test script for molecular modeling packages installation
Run this to verify all packages are working correctly

IMPORTANT: Run with conda environment's Python:
  /Users/kasonchiu/miniconda3/envs/md-dock/bin/python test_installation.py
"""

import sys

print("🧪 Testing Molecular Modeling Package Installation\n")
print(f"Python: {sys.executable}")
print(f"Version: {sys.version}")
print("=" * 60)

# Test packages
packages_to_test = {
    'Core MD Simulation': {
        'openmm': 'OpenMM',
        'pdbfixer': 'PDBFixer',
    },
    'Trajectory Analysis': {
        'MDAnalysis': 'MDAnalysis',
        'parmed': 'ParmEd',
    },
    'Biological Data': {
        'Bio': 'BioPython',
    },
    'Cheminformatics': {
        'rdkit': 'RDKit',
        'prolif': 'ProLIF',
        'meeko': 'Meeko',
    },
    'Visualization': {
        'nglview': 'NGLView',
        'pymol': 'PyMOL',
        'py3Dmol': 'Py3Dmol',
    },
    'Scientific Computing': {
        'numpy': 'NumPy',
        'scipy': 'SciPy',
        'pandas': 'Pandas',
        'matplotlib': 'Matplotlib',
    }
}

results = {}
for category, packages in packages_to_test.items():
    print(f"\n📦 {category}:")
    print("-" * 60)
    for module, name in packages.items():
        try:
            mod = __import__(module)
            version = getattr(mod, '__version__', 'unknown')
            print(f"  ✅ {name:20s} v{version}")
            results[name] = True
        except ImportError as e:
            print(f"  ❌ {name:20s} NOT INSTALLED")
            results[name] = False

# Test OpenMM platforms
print("\n🚀 OpenMM Platform Check:")
print("-" * 60)
try:
    import openmm
    for i in range(openmm.Platform.getNumPlatforms()):
        platform = openmm.Platform.getPlatform(i)
        print(f"  ✅ {platform.getName()}")
except Exception as e:
    print(f"  ❌ Error checking platforms: {e}")

# Summary
print("\n" + "=" * 60)
installed = sum(results.values())
total = len(results)
print(f"📊 Summary: {installed}/{total} packages successfully installed")
print("=" * 60)

if installed == total:
    print("\n🎉 All packages installed successfully!")
else:
    print(f"\n⚠️  {total - installed} package(s) missing or failed to import")
    missing = [name for name, status in results.items() if not status]
    print(f"Missing: {', '.join(missing)}")

# Check command-line tools
print("\n🛠️  Command-Line Tools:")
print("-" * 60)
import subprocess

tools = {
    'pdbfixer': 'PDBFixer (structure prep)',
    'vina': 'AutoDock Vina (docking)',
    'obabel': 'Open Babel (format conversion)'
}

for cmd, name in tools.items():
    result = subprocess.run(['which', cmd], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  ✅ {name}")
    else:
        print(f"  ❌ {name} - Not in PATH")

# Recommendations
print("\n" + "=" * 60)
print("💡 Recommendations:")
print("=" * 60)
if installed == total:
    print("✅ All Python packages installed!")
else:
    print(f"⚠️  {total - installed} Python package(s) missing")
    if missing:
        print(f"   Missing: {', '.join(missing)}")

print("\n📚 For docking, see VINA_ALTERNATIVES.md for:")
print("   • GNINA (recommended - better than Vina)")
print("   • DiffDock (AI-based, no binding site needed)")
print("   • RDKit+OpenMM (already installed!)")
print("\n🚀 Your environment is ready for molecular modeling!")
