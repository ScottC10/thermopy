from thermopy.constants import R
from thermopy.eos.base import EoS
from thermopy.species.species import SPECIES_DATA
import numpy as np




class VanDerWaals(EoS):
    #pure numbers in Van der Waals correlation
    Omega = 1/8
    Psi = 27/64

    def __init__(self, species: str):
        super().__init__(species)
        #Pull required data from species database
        self.Tc = SPECIES_DATA[species]["Tc"]
        self.Pc = SPECIES_DATA[species]["Pc"]

        #coefficients for cubic equation
        self.a = self.Psi   * (R**2 * self.Tc**2)/self.Pc
        self.b = self.Omega * (R*self.Tc)/self.Pc

    #Van der Waals cubic equation in volume: V^3 - (b+RT/P)V^2 + aV/P - ab/P = 0
    def solve(self,T,P):
        #coefficients of cubic polynomial to solve
        coeff = [
            1,
            -(self.b + (R * T) / P),
            self.a / P,
            -((self.a * self.b) / P)
        ]

        """
        Will have three roots:
            -T>Tc  --> 1 real root       : vapour phase molar volume
            -T=Tc  --> 3 identical roots : all represent vapour phase volume
            -T<Tc  --> 3 real roots      : two phase region smaller root is liquid volume,
                                           larger is vapour volume, middle is unphysical
        """
        roots = np.roots(coeff)
        #keep only real roots
        ##print(roots)
        roots = roots[np.abs(roots.imag)<1e-10]
        roots = roots.real

        if len(roots) == 1:
            return {"vapour" : roots[0]}
        elif len(roots)==3:
            roots.sort()
            ##print(roots)

            return {
                "liquid" : roots[0],
                "vapour" : roots[-1]
            }




