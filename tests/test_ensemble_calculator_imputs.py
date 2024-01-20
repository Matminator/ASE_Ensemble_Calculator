import pytest
import warnings

from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT

from ase_ensemble_calculator.ensemble_calculator import Ensemble_Calculator as ES


def test_empty_input():
    with pytest.raises(ValueError):
        ES([])

def test_non_ase_calculator_input():
    with pytest.warns(UserWarning):
        ES([1,2,3])