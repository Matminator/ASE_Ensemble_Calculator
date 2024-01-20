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

def test_1_EMT_setup():

    calcs = [EMT(), EMT(), EMT()] # 3 EMT calculators which are ASE calculator objects
    ensamble = ES(calcs) # Setup Ensemble_Calculator, which sould work for ASE calculators

    assert not isinstance(ensamble, int) # Is not int
    assert isinstance(ensamble, ES) # Is Ensemble_Calculator instance
    assert isinstance(ensamble, Calculator) # Is ASE Calculator instance
    assert ensamble.num_calculators == 3 # Has 3 calculators



    



