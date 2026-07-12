from thermopy.eos.base import EoS, CubicEoS
from thermopy.species.species import get_species
from thermopy.constants import R
import numpy as np

class PengRobinson(CubicEoS):
    OMEGA = 0.07780
    PSI = 0.45724
    sigma = 1 + np.sqrt(2)
    epsilon = 1 - np.sqrt(2)

    def __init__(self, species : str):
        super().__init__(species)
        cur_species = get_species(species)
        self.Tc = cur_species.Tc
        self.Pc = cur_species.Pc
        self.acentric = cur_species.omega

    def P_from_TV(self, T, V) -> float:
        a_alpha, b = self._cubic_parameters(T)
        return (R * T)/(V - b) - a_alpha / (V**2 + 2*b*V - b**2)

    def _cubic_parameters(self, T) -> tuple[float, float]:
        a = self.PSI * (R ** 2 * self.Tc ** 2) / self.Pc
        b = self.OMEGA * (R * self.Tc) / self.Pc

        """peng robinson temperature dependent function alpha(omega, Tr)"""
        k = 0.37464 + 1.54226 * self.acentric - 0.26992 * self.acentric ** 2
        alpha = (1 + k * (1 - np.sqrt(T / self.Tc))) ** 2
        return a*alpha, b

    def _da_dT(self, T):
        tr = T/self.Tc
        j = self.PSI * (R**2 * self.Tc**2)/self.Pc
        k = (0.37464 + 1.54226*self.acentric -0.269928*self.acentric**2)

        alpha_prime = (k *(k* tr**(1/2) - k -1)) / (self.Tc * tr**(1/2))

        return j * alpha_prime


