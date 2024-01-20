
import sys
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT

from ase_ensemble_calculator.ensemble_calculator import Ensemble_Calculator as ES


def test_three_calculatrs_setup():

    calcs = [EMT(), EMT(), EMT()] # 3 EMT calculators which are ASE calculator objects
    ensamble = ES(calcs) # Setup Ensemble_Calculator, which sould work for ASE calculators

    assert isinstance(ensamble, ES) # is Ensemble_Calculator instance
    assert isinstance(ensamble, Calculator)
    assert ensamble.num_calculators == 3



    



