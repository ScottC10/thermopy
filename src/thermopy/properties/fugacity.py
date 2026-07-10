
from thermopy.constants import R
from scipy.integrate import quad
import numpy as np
from thermopy.eos import EoSResult, EoS, CubicEoS
from typing import Union
from thermopy.species import get_species


"""
A generic numerical fugacity coefficient solver using volume explicit departure integral:

    ln(phi) = 1/RT * int(inf->V)[(RT/V)-P]dV + (Z-1) -ln(Z)
    
    Takes an EoSResult object containing all required information along with an equation of state containing
    the correct P_from_RT function
"""


def coefficient(state : EoSResult, EoS : Union[EoS, CubicEoS], cubicEoS=False):
    '''
    A numerical solver for fugacity in the absence of any analytical solution for the used EoS, default solver
    is numerical unless cubicEoS = true
    Uses fugacity coefficient definition from chemical potentials: ln(phi) = 1/RT * int(inf->V)[(RT/V)-P]dV + (Z-1) -ln(Z)
    :param state: EoS result at a T and P
    :param EoS: the EoS used to generate the T and P
    :param cubicEoS: bool deciding which solve method is used default = False
    :return: fugacity coefficient
    '''
    species = get_species(EoS.species)



    if cubicEoS:
        """
        General cubic EoS result taken from: Introduction to Chemical Engineering Thermodynamics 
                                             (J.M. Smith, H.C. Van Ness, M.M. Abbott, M.T. Swihart)
        
        """
        #state parameters
        T = state.T
        P = state.P
        Z = state.Z

        #Cubic EoS parameters
        a, b = EoS._cubic_parameters(T)
        beta = (b * P) / (R * T)
        q = a / (b * R * T)
        I = (1/EoS.sigma - EoS.epsilon) * np.log((Z + EoS.sigma * beta) / (Z + EoS.epsilon * beta))

        ln_phi = Z - 1
        ln_phi += -np.log(Z-beta)
        ln_phi += -q * I

        return np.exp(ln_phi)





    T = state.T
    V = state.V
    Z = state.Z
    P = state.P
    V_ideal = (R*T)/P
    V_inf = 10000*V_ideal

    """Integrates over log space to allow better resolution for liquid phase where contributions matter at small vol"""
    def integrand_log(x):
        V=np.exp(x)
        return float(((R*T/V) - EoS.P_from_TV(T,V))*V)

    integral_log, _ = quad(
        integrand_log,
        np.log(V),
        np.log(V_inf)
    )
    ln_phi = -integral_log/(R * T)
    ln_phi += (Z-1) - np.log(Z)

    #debug prints
    # print("Integral:", -integral_log / (R * T))
    # print("Z-1:", Z - 1)
    # print("-lnZ:", -np.log(Z))
    # print("lnphi:", ln_phi)

    return np.exp(ln_phi)


