import numpy as np

from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.calculator import Calculator

from ase_ensemble_calculator.ensemble_calculator import Ensemble_Calculator as ES

# Mock calculator for testting purpuses
class Noisy_EMT_Calculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, energy = 0, force = 0, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.energy = energy
        self.force = force

    def calculate(self, atoms, properties, system_changes):
        super().calculate(atoms, properties, system_changes)

        energy = self.__calculate_potential_energy(self.atoms)
        forces = self.__calculate_forces(self.atoms)

        atoms.calc = EMT()
        energy += atoms.get_potential_energy()
        forces += atoms.get_forces()

        self.results = {'energy': energy, 'forces': forces}

    def __calculate_potential_energy(self, atoms):
        return self.energy
    
    def __calculate_forces(self, atoms):
        forces = np.ones((len(atoms), 3)) * self.force
        return forces
    

def test_relaxation():
    # H2O molecule
    d = 0.9575
    t = np.pi / 180 * 104.51
    water = Atoms('H2O',
                positions=[(d, 0, 0),
                            (d * np.cos(t), d * np.sin(t), 0),
                            (0, 0, 0)])

    # Ensamble calculator:
    calculators = [EMT(), EMT()]
    ensemble = ES(calculators)
    water.calc = ensemble

    dyn = BFGS(water)
    assert dyn.run(fmax=0.05) # Assert that relaxation has complited

    # Assert that force variances are zero
    assert np.allclose(ensemble.get_forces_standard_deviations(), 0)


def test_relaxation_with_noisy_EMT_calculator():

    # Making 3 tests, as random EMT-calc is added
    for _ in range(3):
        # H2O molecule
        d = 0.9575
        t = np.pi / 180 * 104.51
        water = Atoms('H2O',
                    positions=[(d, 0, 0),
                                (d * np.cos(t), d * np.sin(t), 0),
                                (0, 0, 0)])

        # Ensamble calculator:
        calculators = [EMT(), EMT(), 
                    Noisy_EMT_Calculator(
                            energy = np.random.rand() * .025 - 0.05,
                            force = np.random.rand() * .025 - 0.05)]
        ensemble = ES(calculators)
        water.calc = ensemble

        dyn = BFGS(water)
        assert dyn.run(fmax=0.05) # Assert that relaxation has complited

        # Assert that force variances are diffrent from zero
        assert np.all(ensemble.get_forces_standard_deviations() > 0)


