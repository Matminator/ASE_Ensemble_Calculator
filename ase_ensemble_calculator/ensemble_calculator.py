"""ASE ensemble calculator with ensemble-based uncertainty estimates."""

from __future__ import annotations

from collections.abc import Sequence
import warnings

import numpy as np
from ase.calculators.calculator import Calculator


class Ensemble_Calculator(Calculator):
    """Average energy and forces across multiple ASE calculators.

    The calculator returns the mean energy and forces from the provided
    calculators and can optionally store ensemble-based variance estimates for
    downstream active-learning workflows.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        calculators: Sequence[Calculator],
        compute_variances: bool = True,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        calculators = list(calculators)
        num_calculators = len(calculators)

        if num_calculators == 0:
            raise ValueError("Provided list of calculators is empty.")
        if num_calculators == 1:
            raise ValueError(
                "Provided list of calculators contains only one calculator; "
                "an ensemble calculator is therefore unnecessary."
            )

        non_ase_calculators = sum(not isinstance(calc, Calculator) for calc in calculators)
        if non_ase_calculators > 0:
            warnings.warn(
                f"{non_ase_calculators} out of {len(calculators)} provided calculators "
                "are not ASE calculators. This may cause the ensemble calculator to fail.",
                UserWarning,
            )

        self.calculators = calculators
        self.compute_variances = compute_variances
        self.num_calculators = num_calculators
        self.potential_energy_variance = None
        self.forces_variances = None

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=None):
        """Run all calculators and store mean energy and forces."""
        super().calculate(atoms, properties, system_changes)

        energy = self._calculate_potential_energy(self.atoms)
        forces = self._calculate_forces(self.atoms)
        self.results = {"energy": energy, "forces": forces}

    def get_potential_energy_variance(self):
        """Return the variance of ensemble energies from the latest evaluation."""
        return self.potential_energy_variance

    def get_potential_energy_standard_deviation(self):
        """Return the standard deviation of ensemble energies from the latest evaluation."""
        return np.sqrt(self.potential_energy_variance)

    def get_forces_variances(self):
        """Return per-atom force variances from the latest evaluation."""
        return self.forces_variances

    def get_forces_standard_deviations(self):
        """Return per-atom force standard deviations from the latest evaluation."""
        return np.sqrt(self.forces_variances)

    def _calculate_potential_energy(self, atoms):
        calc_energies = []
        for calc in self.calculators:
            atoms_copy = atoms.copy()
            atoms_copy.calc = calc
            calc_energies.append(atoms_copy.get_potential_energy())

        calc_energies = np.asarray(calc_energies, dtype=float)
        mean_energy = np.mean(calc_energies)

        if self.compute_variances:
            self.potential_energy_variance = np.var(calc_energies)

        return mean_energy

    def _calculate_forces(self, atoms):
        all_forces = []
        for calc in self.calculators:
            atoms_copy = atoms.copy()
            atoms_copy.calc = calc
            all_forces.append(atoms_copy.get_forces())

        all_forces = np.asarray(all_forces, dtype=float)
        mean_forces = np.mean(all_forces, axis=0)

        if self.compute_variances:
            squared_deviation = np.square(all_forces - mean_forces)
            summed_components = np.sum(squared_deviation, axis=0)
            self.forces_variances = np.sum(summed_components, axis=1) / (3 * self.num_calculators)

        return mean_forces
