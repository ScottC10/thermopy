from abc import ABC, abstractmethod

from pygments.unistring import Pc

from thermopy.constants import R
from .results import EoSResult
import numpy as np
from thermopy.species import SPECIES_DATA

class EoS(ABC):
    def __init__(self, species:str):
        '''
            Base class for all equations of state
            :param species: species to be analysed
        '''
        self.species = species

    def _build_result(self, T, P, Z, phase):
        '''
        constructs a dataclass containing the thermodynamic information at a specific T and P
        :param T: temperature [K]
        :param P: pressure [Pa]
        :param Z: compressibility factor
        :param phase: phase of the species
        :return:
        '''
        V = (Z* R * T) / P
        return EoSResult(
            T=T, P=P, V=V, Z=Z,
            phase=phase,
            eos_name=type(self).__name__
        )


    """
    Abstract Methods 
    """
    @abstractmethod
    def solve(self, T, P) ->list[EoSResult] :
        '''

        :param T: temperature [K]
        :param P: pressure [P]
        :return: compressibility
        '''
        pass
    @abstractmethod
    def P_from_TV(self, T, V) -> float:
        '''
        finds the pressure of the species at a given temperature and pressure
        :param T: temperature [K]
        :param V: volume [m^3]
        :return: pressure [Pa]
        '''
        pass


class CubicEoS(EoS):
    PSI = None
    OMEGA = None

    sigma = None
    epsilon = None
    def solve(self, T, P):
        '''
        Solves cubic EoS at specific T and P returns a number of roots dependent on phase
        :param T: temperature [K]
        :param P: pressure [Pa]
        :return: compressibility factor roots of cubic EoS
        '''
        """
         Will have three roots:
            -T>Tc  --> 1 real root       : vapour phase molar volume
            -T=Tc  --> 3 identical roots : all represent vapour phase volume
            -T<Tc  --> 3 real roots      : two phase region smaller root is liquid volume,
                                                   larger is vapour volume, middle is unphysical
        """
        a, b = self._cubic_parameters(T)

        A = (a * P) / (R ** 2 * T ** 2)  # Dimensionless quantities
        B = (b * P) / (R * T)

        coeffs = [
            1,
            ((self.epsilon + self.sigma - 1) * B - 1),
            (A - (self.epsilon + self.sigma) * B + (self.epsilon * self.sigma - self.epsilon - self.sigma) * B ** 2),
            -(A * B + self.epsilon * self.sigma * B ** 2 * (1 + B))
        ]
        roots = np.roots(coeffs)
        roots = roots[np.abs(roots.imag) < 1e-10]
        roots = roots.real
        roots.sort()
        if len(roots) not in (1, 3):
            raise ValueError(f"Unexpected number of real rots: {len(roots)}")
        return [self._build_result(T, P, Z, phase=None) for Z in roots]



    @abstractmethod
    def _cubic_parameters(self, T) -> tuple[float, float]:
        """
        finds the cubic parameters required for the cubic equations of state
        :param T: temperature [T]
        :return: (a, b) or (a*alpha, b) in the case of temperature dependent EoS at given T
        """


