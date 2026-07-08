from abc import ABC, abstractmethod
from thermopy.constants import R
from .results import EoSResult

class EoS(ABC):
    def __init__(self, species:str):
        '''
            Base class for all equations of state
            :param species: species to be analysed
        '''

        self.species = species


    def build_result(self, T, P, Z, phase):
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
    @abstractmethod
    def cubic_parameters(self, T) -> tuple[float, float]:
        """
        finds the cubic parameters required for the cubic equations of state
        :param T: temperature [T]
        :return: (a, b) or (a*alpha, b) in the case of temperature dependent EoS at given T
        """
