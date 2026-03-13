# ASE Ensemble Calculator

[![Tests](https://github.com/Matminator/ASE_Ensemble_Calculator/actions/workflows/pytests.yaml/badge.svg)](https://github.com/Matminator/ASE_Ensemble_Calculator/actions/workflows/pytests.yaml)

ASE calculator wrapper for ensemble-based uncertainty estimates in atomistic simulation and active-learning workflows.

## Why this exists

When working with ensembles of interatomic models, it is often useful to:
- average energies and forces across several ASE calculators
- expose uncertainty estimates from disagreement within the ensemble
- plug the result directly into existing ASE-based simulation or relaxation workflows

This package provides a lightweight calculator that does exactly that.

## Key features

- wraps multiple ASE calculators in a single calculator interface
- returns mean energy and mean forces across the ensemble
- stores ensemble variance estimates for energy and per-atom forces
- integrates directly with standard ASE workflows, including structure relaxation
- tested with pytest across Windows, macOS, and Linux in GitHub Actions

## Installation

Install from a local checkout:

```bash
pip install .
```

Install with development dependencies:

```bash
pip install -r requirements_dev.txt
pip install -e .
```

## Basic usage

```python
from ase_ensemble_calculator import EnsembleCalculator

calc1 = calculator1()  # ASE calculator
calc2 = calculator2()  # Another ASE calculator

ensemble = EnsembleCalculator([calc1, calc2])
atoms.calc = ensemble

energy = atoms.get_potential_energy()
forces = atoms.get_forces()

energy_variance = ensemble.get_potential_energy_variance()
forces_variances = ensemble.get_forces_variances()
```

## API summary

The main public class is `EnsembleCalculator`.

It:
- takes a sequence of ASE calculators
- returns mean ensemble predictions for `energy` and `forces`
- optionally computes variance estimates via `compute_variances=True`

Available uncertainty helpers:
- `get_potential_energy_variance()`
- `get_potential_energy_standard_deviation()`
- `get_forces_variances()`
- `get_forces_standard_deviations()`

## Testing

Run the test suite with:

```bash
pytest tests/
```

The repository includes tests for:
- basic ensemble behavior on ASE EMT calculators
- deterministic behavior for simple mock calculators
- variance and standard-deviation outputs
- compatibility with ASE relaxation workflows

## Repository structure

- `ase_ensemble_calculator/ensemble_calculator.py`: calculator implementation
- `tests/`: pytest-based test suite
- `.github/workflows/pytests.yaml`: CI configuration
- `ASE_Ensemble_Calculator.ipynb`: exploratory notebook / usage notes

## Intended use

This package was developed for ensemble-based active-learning workflows in atomistic simulation, but the calculator is reusable anywhere multiple ASE-compatible calculators should be aggregated into a single interface.
