#!/usr/bin/env python3
"""
Minimal OpenMM example.

What it does:
- Creates a short 2-residue peptide PDB programmatically
- Builds a system using an Amber force field and runs a very short (few ps) integrator
- Writes a tiny trajectory (`traj.dcd`) and `log.txt` if OpenMM is installed

Note: This is illustrative. For production runs, prepare full structures, solvate, add ions, equilibrate, and run longer.
"""
import os
import sys

def run_openmm_test():
    try:
        from openmm import app
        import openmm as mm
        from openmm import unit
    except Exception:
        print('OpenMM not available. Install OpenMM (conda recommended).')
        return

    # Minimal Lennard-Jones-like test (no force field needed). This creates
    # a small cubic lattice of neutral particles interacting via a
    # NonbondedForce (LJ) and runs a short Langevin dynamics to validate OpenMM.
    system = mm.System()
    n_per_side = 3
    n_particles = n_per_side**3
    mass = 39.948 * unit.amu  # approximate argon mass
    for i in range(n_particles):
        system.addParticle(mass)

    nonbond = mm.NonbondedForce()
    nonbond.setNonbondedMethod(mm.NonbondedForce.PME)  # use PME/periodic
    sigma = 0.34 * unit.nanometer
    epsilon = 0.238 * unit.kilojoule_per_mole
    for i in range(n_particles):
        nonbond.addParticle(0.0, sigma, epsilon)
    system.addForce(nonbond)

    # Box vectors: 2.0 nm cube
    box_size = 2.0 * unit.nanometer
    system.setDefaultPeriodicBoxVectors((box_size, 0, 0), (0, box_size, 0), (0, 0, box_size))

    integrator = mm.LangevinIntegrator(120*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    platform = mm.Platform.getPlatformByName('CPU')
    simulation = app.Simulation(app.Topology(), system, integrator, platform)

    # Create cubic lattice positions
    positions = []
    spacing = box_size / (n_per_side + 1)
    for x in range(1, n_per_side+1):
        for y in range(1, n_per_side+1):
            for z in range(1, n_per_side+1):
                positions.append(mm.Vec3(x*spacing, y*spacing, z*spacing))
    simulation.context.setPositions(positions)

    simulation.reporters.append(app.StateDataReporter('log.txt', 100, step=True, potentialEnergy=True, temperature=True))
    print('Running a short OpenMM test (500 steps) ...')
    simulation.step(500)
    print('Done. Check `log.txt` for energy/temperature information.')

if __name__ == '__main__':
    run_openmm_test()
