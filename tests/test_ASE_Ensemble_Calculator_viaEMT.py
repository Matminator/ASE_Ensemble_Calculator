
from ase.calculators.emt import EMT



def test_three_calculatrs_setup():

    calcs = [EMT(), EMT(), EMT()]

    from ASE_Ensemble_Calculator import ASE_Ensemble_Calculator as EC

    ensamble = EC(calcs)


    



