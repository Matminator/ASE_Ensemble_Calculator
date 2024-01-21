# ASE Ensemble Calculator

Simple and streamlined ASE calculator for machine learning ensemble and active learning methods.

Basic usage:
```
    calc1 = calculator1() # A ASE calculator
    calc2 = calculator2() # A diffrent ASE calculator

    ensamble = ES([calc1, calc2])
    atoms.calc = ensamble

    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    energy_variance = ensamble.get_potential_energy_variance()
    forces_variances = ensamble.get_forces_variances()
```
