from pip._internal.models import selection_prefs

from thermopy.constants import R
from thermopy.eos.base import EoS
from thermopy.species.species import SPECIES_DATA
import numpy as np




class VanDerWaals(EoS):
    #pure numbers in Van der Waals correlation
    OMEGA = 1/8
    PSI = 27/64

    sigma = 0
    epsilon = 0

    def __init__(self, species: str):
        super().__init__(species)
        #Pull required species data from the species database
        self.Tc = SPECIES_DATA[species]["Tc"]
        self.Pc = SPECIES_DATA[species]["Pc"]


    def solve(self,T,P):
        """
        Will have three roots:
            -T>Tc  --> 1 real root       : vapour phase molar volume
            -T=Tc  --> 3 identical roots : all represent vapour phase volume
            -T<Tc  --> 3 real roots      : two phase region smaller root is liquid volume,
                                           larger is vapour volume, middle is unphysical
        """
        a, b = self.cubic_parameters(T)

        A = (a * P) / (R**2 * P**2) #Dimensionless quantities
        B = (b * P) / (R * T)

        coeffs = [
            1,
            ((self.epsilon + self.sigma - 1) * B -1),
            (A - (self.epsilon+self.sigma)*B + (self.epsilon*self.sigma -self.epsilon-self.sigma)*B**2),
            (A*B + self.epsilon*self.sigma*B**2*(1+B))
        ]
        roots = np.roots(coeffs)
        roots = roots[np.abs(roots.imag) < 1e-10]
        roots = roots.real
        roots.sort()
        if len(roots) not in (1, 3):
            raise ValueError(f"Unexpected number of real rots: {len(roots)}")
        return [self.build_result(T, P, Z, phase=None) for Z in roots]


    def cubic_parameters(self, T) -> tuple[float, float]:
        a = self.PSI * (R**2 * self.Tc**2) / self.Pc
        b = self.OMEGA * (R * self.Tc)/self.Pc
        return a, b

    def P_from_TV(self, T, V) -> float:
        a, b = self.cubic_parameters(T)
        return (R*T)/(V-b) - a/((V+self.epsilon*b)*(V+self.sigma*b))






