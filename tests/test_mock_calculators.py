import numpy as np

from ase import Atoms
from ase.calculators.calculator import Calculator

from ase_ensemble_calculator import EnsembleCalculator as ES


class MockCalculator(Calculator):
    """Simple mock calculator for ensemble calculator tests."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, energy=0, force=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.energy = energy
        self.force = force

    def calculate(self, atoms, properties, system_changes):
        super().calculate(atoms, properties, system_changes)
        energy = self._calculate_potential_energy(self.atoms)
        forces = self._calculate_forces(self.atoms)
        self.results = {"energy": energy, "forces": forces}

    def _calculate_potential_energy(self, atoms):
        return self.energy

    def _calculate_forces(self, atoms):
        return np.ones((len(atoms), 3)) * self.force


def test_simple_calculators_output_shapes():
    calc1 = MockCalculator()
    calc2 = MockCalculator()

    ensemble = ES([calc1, calc2])

    atoms = Atoms(["H", "H"], positions=([0, 0, 0], [0.5, 0.5, 0.5]))
    atoms.calc = ensemble

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    assert energy.shape == ()
    assert forces.shape == (2, 3)


def test_compute_variances_option():
    calc1 = MockCalculator()
    calc2 = MockCalculator()

    ensemble = ES([calc1, calc2], compute_variances=False)

    atoms = Atoms(["H", "H"], positions=([0, 0, 0], [0.5, 0.5, 0.5]))
    atoms.calc = ensemble

    energy1 = atoms.get_potential_energy()
    forces1 = atoms.get_forces()

    assert energy1.shape == ()
    assert forces1.shape == (2, 3)
    assert ensemble.get_potential_energy_variance() is None
    assert ensemble.get_forces_variances() is None

    ensemble = ES([calc1, calc2])
    atoms.calc = ensemble

    energy2 = atoms.get_potential_energy()
    forces2 = atoms.get_forces()

    assert np.allclose(energy1, energy2)
    assert np.allclose(forces1, forces2)
    assert ensemble.get_potential_energy_variance() is not None
    assert not ensemble.get_forces_variances().any(None)
    assert ensemble.get_potential_energy_standard_deviation() is not None
    assert not ensemble.get_forces_standard_deviations().any(None)


def test_simple_calculators_1():
    calc1 = MockCalculator(energy=-3, force=-3)
    calc2 = MockCalculator(energy=3, force=3)

    ensemble = ES([calc1, calc2])

    atoms = Atoms(["H", "H"], positions=([0, 0, 0], [0.5, 0.5, 0.5]))
    atoms.calc = ensemble

    assert np.allclose(atoms.get_potential_energy(), 0)
    assert np.allclose(atoms.get_forces(), 0)
    assert np.allclose(ensemble.get_potential_energy_variance(), 9)
    assert np.allclose(ensemble.get_forces_variances(), 9)
    assert np.allclose(ensemble.get_potential_energy_standard_deviation(), 3)
    assert np.allclose(ensemble.get_forces_standard_deviations(), 3)


def test_simple_calculators_2():
    pos = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float64)
    atoms = Atoms("H2C3", positions=pos)

    for _ in range(5):
        num_calculators = np.random.randint(2, 10)
        calcs = []
        for _ in range(num_calculators):
            calc = MockCalculator(energy=np.random.rand() * 10 - 5, force=np.random.rand() * 10 - 5)
            calcs.append(calc)

        all_energies = []
        average_energy = 0
        for calc in calcs:
            atoms.calc = calc
            energy = atoms.get_potential_energy()
            all_energies.append(energy)
            average_energy += energy
        average_energy /= num_calculators

        all_forces = []
        average_forces = np.zeros((len(atoms), 3))
        for calc in calcs:
            atoms.calc = calc
            forces = atoms.get_forces()
            all_forces.append(forces)
            average_forces += forces
        average_forces /= num_calculators

        ensemble = ES(calcs)
        atoms.calc = ensemble
        ensemble_energy = atoms.get_potential_energy()
        ensemble_forces = atoms.get_forces()

        assert np.allclose(ensemble_energy, average_energy)
        assert np.allclose(ensemble_forces, average_forces)
        assert np.allclose(ensemble.get_potential_energy_variance(), np.var(all_energies))
        assert np.allclose(ensemble.get_forces_variances(), np.var(all_forces))
