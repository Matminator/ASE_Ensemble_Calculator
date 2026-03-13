"""Minimal example for using the ASE ensemble calculator with ASE EMT."""

from ase import Atoms
from ase.calculators.emt import EMT

from ase_ensemble_calculator import EnsembleCalculator


def main():
    atoms = Atoms("H2O", positions=[(0.9575, 0.0, 0.0), (-0.2390, 0.9270, 0.0), (0.0, 0.0, 0.0)])

    ensemble = EnsembleCalculator([EMT(), EMT(), EMT()])
    atoms.calc = ensemble

    print(f"Mean energy: {atoms.get_potential_energy():.6f} eV")
    print(f"Mean forces:\n{atoms.get_forces()}")
    print(f"Energy variance: {ensemble.get_potential_energy_variance():.6e}")
    print(f"Force variances: {ensemble.get_forces_variances()}")


if __name__ == "__main__":
    main()
