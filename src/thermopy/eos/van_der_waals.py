from thermopy.constants import R
from thermopy.eos.base import EoS, CubicEoS
from thermopy.species import get_species
import numpy as np




class VanDerWaals(CubicEoS):
    #pure numbers in Van der Waals correlation
    OMEGA = 1/8
    PSI = 27/64

    sigma = 0
    epsilon = 0

    def __init__(self, species: str):
        super().__init__(species)
        #Pull required species data from the species database
        cur_species = get_species(species)
        self.Tc = cur_species.Tc
        self.Pc = cur_species.Pc

    def _cubic_parameters(self, T) -> tuple[float, float]:
        a = self.PSI * (R**2 * self.Tc**2) / self.Pc
        b = self.OMEGA * (R * self.Tc)/self.Pc
        return a, b

    def P_from_TV(self, T, V) -> float:
        a, b = self._cubic_parameters(T)
        return (R*T)/(V-b) - a/((V+self.epsilon*b)*(V+self.sigma*b))

    def _da_dT(self, T):
        #VdW does not have a temperature dependence
        pass





