import numpy as np
from ase.calculators.calculator import Calculator
import warnings

class Ensemble_Calculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, calculators: list, compute_variances = True, *args, **kwargs):
        super().__init__(*args, **kwargs)

        num_calculators = len(calculators)

        # Testing input:
        if num_calculators == 0:
            raise ValueError('Provided list of calculators is empty (length is 0)')
        
        # Testing type of list
        nun_ASE_calcs = 0
        for calc in calculators:
            if not isinstance(calc, Calculator):
                nun_ASE_calcs += 1
        if nun_ASE_calcs > 0:
            warnings.warn(
                f"{nun_ASE_calcs} out of {len(calculators)} elements of the provided calculators-list "
                "are not ASE calculators. This may result in the failure of this calculator.",
                UserWarning
            )

        self.calculators = calculators
        self.potential_energy_variance = None
        self.forces_variances = None
        self.num_calculators = num_calculators
        self.compute_variances = compute_variances

    def calculate(self, atoms, properties, system_changes):
        super().calculate(atoms, properties, system_changes)

        energy = self.__calculate_potential_energy(self.atoms)
        forces = self.__calculate_forces(self.atoms)
        self.results = {'energy': energy, 'forces': forces}

    def get_potential_energy_variance(self):
        return self.potential_energy_variance
    
    def get_potential_energy_standard_deviation(self):
        return np.sqrt(self.potential_energy_variance)
    
    def get_forces_variances(self):
        return self.forces_variances

    def get_forces_standard_deviations(self):
        return np.sqrt(self.forces_variances)


    # Private methods for computation of energy and forces:

    def __calculate_potential_energy(self, atoms):

        calc_energies = []
        for calc in self.calculators:
            atoms_copy = atoms.copy()
            atoms_copy.calc = calc
            calc_energies.append(atoms_copy.get_potential_energy())
        
        mean_energy = np.mean(calc_energies)
    
        if  self.compute_variances:    
            # Computing the variance of the energies:
            self.potential_energy_variance = np.var(calc_energies)

        return mean_energy
    
    def __calculate_forces(self, atoms):

        all_forces = []
        for calc in self.calculators:
            atoms_copy = atoms.copy()
            atoms_copy.calc = calc
            all_forces.append(atoms_copy.get_forces())
        
        mean_forces = np.mean(all_forces, axis = 0)
        
        if  self.compute_variances:
            # Computing the variance of all forces:
            x = (all_forces - mean_forces)**2

            x = np.sum(x, axis = 0)
            x = np.sum(x, axis = 1) / (3 * self.num_calculators)
            self.forces_variances = np.sqrt(x) 

        return mean_forces