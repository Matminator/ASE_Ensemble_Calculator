import numpy as np
from ase.calculators.calculator import Calculator
import warnings

class  ASE_Ensemble_Calculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, calculators: list):

        num_models = len(calculators)

        # Testing input:
        if num_models == 0:
            raise ValueError('Provided list of calculators is empty (length is 0)')
        
        # Testing type of list
        nun_ASE_calcs = 0
        for calc in calculators:
            if type(calc) != Calculator:
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
        self.num_models = num_models

    def calculate(self, atoms, properties, system_changes):
        # Check which properties need to be calculated
        energy = 0.0
        forces = None

        if 'energy' in properties:
            # Perform energy calculation here (replace this with your actual calculation)
            energy = self.__calculate_potential_energy(atoms)

        if 'forces' in properties:
            # Perform forces calculation here (replace this with your actual calculation)
            forces = self.__calculate_forces(atoms)

        # Store the calculated values
        self.results = {'energy': energy, 'forces': forces}

    def get_potential_energy_variance(self):
        return self.potential_energy_variance
    
    def get_potential_energy_standard_deviation(self):
        return np.sqrt(self.potential_energy_variance)
    
    def get_forces_variances(self):
        pass

    def get_forces_standard_deviations(self):
        pass

    def __calculate_potential_energy(self, atoms):

        calc_energies = []
        for calc in self.calculators:
            atoms_copy = atoms.copy()
            atoms.calc = calc
            calc_energies.append(atoms.get_potential_energy())
        

        
        return 0
    
    def __calculate_forces(self, atoms):

        calc_forces = []
        for calc in self.calculators:
            atoms_copy = atoms.copy()
            atoms.calc = calc
            calc_forces.append(atoms.get_forces())
        
        

        return 0

