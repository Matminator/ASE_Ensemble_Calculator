import numpy as np

from ase import Atoms
from ase.calculators.calculator import Calculator

from ase_ensemble_calculator.ensemble_calculator import Ensemble_Calculator as ES


# Mock calculator for testting purpuses
class mock_calculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, energy = 0, force = 0, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.energy = energy
        self.force = force

    def calculate(self, atoms, properties, system_changes):
        super().calculate(atoms, properties, system_changes)

        energy = self.__calculate_potential_energy(self.atoms)
        forces = self.__calculate_forces(self.atoms)
        self.results = {'energy': energy, 'forces': forces}

    def __calculate_potential_energy(self, atoms):
        return self.energy
    
    def __calculate_forces(self, atoms):
        forces = np.ones((len(atoms), 3)) * self.force
        return forces



def test_simple_calculators_output_shapes():

    calc1 = mock_calculator()
    calc2 = mock_calculator()

    ensamble = ES([calc1, calc2])

    atoms = Atoms(["H", "H"], positions= ([0, 0, 0], [0.5, 0.5, 0.5]))
    atoms.calc = ensamble

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    assert energy.shape == ()
    assert forces.shape == (2, 3)

def test_compute_variances_option():

    calc1 = mock_calculator()
    calc2 = mock_calculator()

    ensamble = ES([calc1, calc2], compute_variances = False)

    atoms = Atoms(["H", "H"], positions= ([0, 0, 0], [0.5, 0.5, 0.5]))
    atoms.calc = ensamble

    energy1 = atoms.get_potential_energy()
    forces1 = atoms.get_forces()

    # Testing what energy and forces have been computed:
    assert energy1.shape == ()
    assert forces1.shape == (2, 3)

    # Testing that variances have not been computed:
    assert ensamble.get_potential_energy_variance() == None
    assert ensamble.get_forces_variances() == None

    # Ensable with variance computation:
    ensamble = ES([calc1, calc2])
    atoms.calc = ensamble

    energy2 = atoms.get_potential_energy()
    forces2 = atoms.get_forces()

    # Testing that the variance computations does not effect the energy and forces:
    assert np.allclose(energy1, energy2)
    assert np.allclose(forces1, forces2)

    # Testing that variances have been computed:
    assert not ensamble.get_potential_energy_variance() == None
    assert not ensamble.get_forces_variances().any(None)
    # Testing that standard divations may be computed:
    assert not ensamble.get_potential_energy_standard_deviation() == None
    assert not ensamble.get_forces_standard_deviations().any(None)

def test_simple_calculators_1():

    calc1 = mock_calculator(energy = -3, force = -3)
    calc2 = mock_calculator(energy = 3, force = 3)

    ensamble = ES([calc1, calc2])

    atoms = Atoms(["H", "H"], positions= ([0, 0, 0], [0.5, 0.5, 0.5]))
    atoms.calc = ensamble

    # Testing resulting energy and forces:
    assert np.allclose(atoms.get_potential_energy(), 0)
    assert np.allclose(atoms.get_forces(), 0)       

    # Testing variances:
    assert np.allclose(ensamble.get_potential_energy_variance(), 9)
    assert np.allclose(ensamble.get_forces_variances(), 9)

    # Testing standard diviations:
    assert np.allclose(ensamble.get_potential_energy_standard_deviation(), 3)
    assert np.allclose(ensamble.get_forces_standard_deviations(), 3)

def test_simple_calculators_2():

    pos = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=np.float64)
    atoms = Atoms('H2C3', positions = pos)

    # Doing tests N times:
    N = 5
    for _ in range(N):

        # Creating M calculators:
        M = np.random.randint(2, 10)
        calcs = []
        for _ in range(M):
            calc = mock_calculator(energy = np.random.rand() * 10 - 5, force = np.random.rand() * 10 - 5)
            calcs.append(calc)

        # Average energy:
        all_energies = []
        avg_energy = 0
        for calc in calcs:
            atoms.calc = calc
            energy = atoms.get_potential_energy()
            all_energies.append(energy)
            avg_energy += energy
        avg_energy /= M

        # Average forces:
        all_forces = []
        avg_forces = np.zeros((len(atoms), 3))
        for calc in calcs:
            atoms.calc = calc
            froces = atoms.get_forces()
            all_forces.append(froces)
            avg_forces += froces
        avg_forces /= M

        # Creating ensemble:
        ensamble = ES(calcs)
        atoms.calc = ensamble
        ens_energy = atoms.get_potential_energy()
        ens_forces = atoms.get_forces()

        # Testing resulting energy and forces:
        assert np.allclose(ens_energy, avg_energy)
        assert np.allclose(ens_forces, avg_forces)

        # Testing variances:
        assert np.allclose(ensamble.get_potential_energy_variance(), np.var(all_energies))
        assert np.allclose(ensamble.get_forces_variances(), np.var(all_forces))

    