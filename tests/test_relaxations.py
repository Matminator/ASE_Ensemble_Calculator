import numpy as np

from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.calculator import Calculator

from ase_ensemble_calculator import EnsembleCalculator as ES


class NoisyEMTCalculator(Calculator):
    """Mock calculator that perturbs EMT energy and forces."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, energy=0, force=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.energy = energy
        self.force = force

    def calculate(self, atoms, properties, system_changes):
        super().calculate(atoms, properties, system_changes)

        energy = self._calculate_potential_energy(self.atoms)
        forces = self._calculate_forces(self.atoms)

        atoms.calc = EMT()
        energy += atoms.get_potential_energy()
        forces += atoms.get_forces()

        self.results = {"energy": energy, "forces": forces}

    def _calculate_potential_energy(self, atoms):
        return self.energy

    def _calculate_forces(self, atoms):
        return np.ones((len(atoms), 3)) * self.force


def test_relaxation():
    d = 0.9575
    t = np.pi / 180 * 104.51
    water = Atoms(
        "H2O",
        positions=[(d, 0, 0), (d * np.cos(t), d * np.sin(t), 0), (0, 0, 0)],
    )

    calculators = [EMT(), EMT()]
    ensemble = ES(calculators)
    water.calc = ensemble

    dyn = BFGS(water)
    assert dyn.run(fmax=0.05)
    assert np.allclose(ensemble.get_forces_standard_deviations(), 0)


def test_relaxation_with_noisy_emt_calculator():
    for _ in range(3):
        d = 0.9575
        t = np.pi / 180 * 104.51
        water = Atoms(
            "H2O",
            positions=[(d, 0, 0), (d * np.cos(t), d * np.sin(t), 0), (0, 0, 0)],
        )

        calculators = [
            EMT(),
            EMT(),
            NoisyEMTCalculator(
                energy=np.random.rand() * 0.025 - 0.05,
                force=np.random.rand() * 0.025 - 0.05,
            ),
        ]
        ensemble = ES(calculators)
        water.calc = ensemble

        dyn = BFGS(water)
        assert dyn.run(fmax=0.05)
        assert np.all(ensemble.get_forces_standard_deviations() > 0)
