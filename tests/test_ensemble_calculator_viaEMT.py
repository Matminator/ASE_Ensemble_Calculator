import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT

from ase_ensemble_calculator import EnsembleCalculator as ES


def test_base_emt_setup():
    calcs = [EMT(), EMT(), EMT()]
    ensemble = ES(calcs)

    assert not isinstance(ensemble, int)
    assert isinstance(ensemble, ES)
    assert isinstance(ensemble, Calculator)
    assert ensemble.num_calculators == 3
    assert ensemble.get_potential_energy_variance() is None
    assert ensemble.get_forces_variances() is None


def test_same_emt_calculations():
    # Test that the ensemble matches EMT when all calculators are identical.
    for _ in range(3):
        pos = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0]], dtype=np.float64)
        pos += np.random.rand(5, 3) * 0.2
        atoms = Atoms("H2C3", positions=pos)

        calcs = [EMT(), EMT(), EMT()]
        ensemble = ES(calcs)

        atoms.calc = ensemble
        ensemble_energy = atoms.get_potential_energy()
        ensemble_forces = atoms.get_forces()

        atoms.calc = EMT()
        emt_energy = atoms.get_potential_energy()
        emt_forces = atoms.get_forces()

        assert np.allclose(ensemble_energy, emt_energy)
        assert np.allclose(ensemble_forces, emt_forces)
        assert np.allclose(ensemble.get_potential_energy_variance(), 0)
        assert np.allclose(ensemble.get_forces_variances(), 0)
        assert np.allclose(ensemble.get_potential_energy_standard_deviation(), 0, atol=1e-4)
        assert np.allclose(ensemble.get_forces_standard_deviations(), 0, atol=1e-4)
