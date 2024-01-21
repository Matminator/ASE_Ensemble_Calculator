import pytest

from ase.calculators.emt import EMT

from ase_ensemble_calculator.ensemble_calculator import Ensemble_Calculator as ES


def test_empty_input():
    with pytest.raises(ValueError):
        ES([])

def test_one_calculator_input():
    with pytest.raises(ValueError):
        ES([EMT()])

def test_non_ase_calculator_input1():
    with pytest.warns(UserWarning):
        ES([1,2,3])

def test_non_ase_calculator_input2():
    with pytest.warns(UserWarning):
        ES([EMT(),EMT(), [], EMT()])