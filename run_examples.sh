#!/bin/bash
# Run MD/Docking examples in the md-dock environment

echo "=== Running OpenMM Example ==="
$HOME/miniconda3/bin/conda run -n md-dock python openmm_example.py
echo ""

echo "=== Running Docking Example ==="
$HOME/miniconda3/bin/conda run -n md-dock python docking_example.py
echo ""

echo "✅ Both examples completed!"
echo ""
echo "Check outputs:"
echo "  - log.txt (OpenMM simulation log)"
echo "  - out.pdbqt (docking poses)"
