from thermopy.eos.base import EoS
from thermopy.species.species import SPECIES_DATA
from thermopy.constants import R


import numpy as np

class PengRobinson(EoS):
    OMEGA = 0.07780
    PSI = 0.45724


    def __init__(self, species : str):
        super().__init__(species)

        self.Tc = SPECIES_DATA[species]["Tc"]
        self.Pc = SPECIES_DATA[species]["Pc"]
        self.acentric = SPECIES_DATA[species]["omega"]



    def solve(self, T, P):
        '''
        solves cubic polynomial in Z
        :param P: pressure [Pa]
        :param T: temperature [K]
        :return: list of Z, size depends on how many phases present
        '''
        a_alpha, b = self.cubic_parameters(T)

        """dimensionless parameters"""
        A = (a_alpha * P) / (R**2 * T**2)
        B = (b * P) / (R * T)

        """cubic coefficients"""
        coeffs = [
            1,
            -(1-B),
            (A - 3*B**2 - 2*B),
            -(A*B - B**2 -B**3)
        ]
        roots = np.roots(coeffs)
        roots = roots[np.abs(roots.imag)<1e-10]
        roots = roots.real
        roots.sort()
        if len(roots) not in (1, 3):
            raise ValueError(f"Unexpected number of real rots: {len(roots)}")
        return [self.build_result(T, P, Z, phase=None) for Z in roots]

    def P_from_TV(self, T, V) -> float:
        a_alpha, b = self.cubic_parameters(T)
        return (R * T)/(V - b) - a_alpha / (V**2 + 2*b*V - b**2)

    def cubic_parameters(self, T) -> tuple[float, float]:
        a = self.PSI * (R ** 2 * self.Tc ** 2) / self.Pc
        b = self.OMEGA * (R * self.Tc) / self.Pc

        """peng robinson temperature dependent function alpha(omega, Tr)"""
        k = 0.37464 + 1.54226 * self.acentric - 0.26992 * self.acentric ** 2
        alpha = (1 + k * (1 - np.sqrt(T / self.Tc))) ** 2
        return a*alpha, b