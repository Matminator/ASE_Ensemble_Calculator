import numpy as np
from ase import Atoms

from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT

from ase_ensemble_calculator.ensemble_calculator import Ensemble_Calculator as ES

def test_base_EMT_setup():

    calcs = [EMT(), EMT(), EMT()] # 3 EMT calculators which are ASE calculator objects
    ensamble = ES(calcs) # Setup Ensemble_Calculator, which sould work for ASE calculators

    assert not isinstance(ensamble, int) # Is not int
    assert isinstance(ensamble, ES) # Is Ensemble_Calculator instance
    assert isinstance(ensamble, Calculator) # Is ASE Calculator instance
    assert ensamble.num_calculators == 3 # Has 3 calculators
    assert ensamble.get_potential_energy_variance() == None
    assert ensamble.get_forces_variances() == None

def test_same_EMT_calculations():

    # Testing if calculations yilds the same results for 3 diffrent positions:
    for _ in range(3):
        pos = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0]], dtype=np.float64)
        pos += np.random.rand(5,3) * 0.2
        atoms = Atoms('H2C3', positions = pos)

        calcs = [EMT(), EMT(), EMT()] # 3 EMT calculators which are ASE calculator objects
        ensamble = ES(calcs) # Setup Ensemble_Calculator, which sould work for ASE calculators

        atoms.calc = ensamble

        ens_energy = atoms.get_potential_energy()
        ens_forces = atoms.get_forces()
        
        atoms.calc = EMT()

        emt_energy = atoms.get_potential_energy()
        emt_forces = atoms.get_forces()

        # Testing if calculations yilds the same results:
        assert np.allclose(ens_energy, emt_energy)
        assert np.allclose(ens_forces, emt_forces)

        # Testing if variance is zero:
        assert np.allclose(ensamble.get_potential_energy_variance(), 0)
        assert np.allclose(ensamble.get_forces_variances(), 0)

        # Testing if standard divations are zero:
        assert np.allclose(ensamble.get_potential_energy_standard_deviation(), 0, atol=1e-4)
        assert np.allclose(ensamble.get_forces_standard_deviations(), 0, atol=1e-4)
    



    



