
import sys
sys.path.append('..')
sys.path.append('/home/matnis/ASE_Ensemble_Calculator/')

from ASE_Ensemble_Calculator.ase_ensemble_calculator.ase_ensemble_calculator import Ensemble_Calculator as ES
from ase.calculators.emt import EMT


def test_three_calculatrs_setup():

    calcs = [EMT(), EMT(), EMT()]

    ensamble = ES(calcs)



    



