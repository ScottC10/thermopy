
from thermopy.constants import R
from scipy.integrate import quad
import numpy as np
from thermopy.eos import EoSResult, CubicEoS

"""
A generic numerical fugacity coefficient solver using volume explicit departure integral:

    ln(phi) = 1/RT * int(inf->V)[(RT/V)-P]dV + (Z-1) -ln(Z)
    
    Takes an EoSResult object containing all required information along with an equation of state containing
    the correct P_from_RT function
"""


def coefficient(state : EoSResult, EoS : CubicEoS ):
    '''
    Uses fugacity coefficient definition from chemical potentials: ln(phi) = 1/RT * int(inf->V)[(RT/V)-P]dV + (Z-1) -ln(Z)
    :param state: EoS result at a T and P
    :param EoS: the EoS used to generate the T and P
    :return: fugacity coefficient
    '''
    T = state.T
    V = state.V
    Z = state.Z
    P = state.P
    V_ideal = (R*T)/P
    V_inf = 10000*V_ideal

    # def integrand(V_prime):
    #     return(
    #             ((R*T)/V_prime) - EoS.P_from_TV(T, V_prime)
    #     )

    # integral, _ = quad(integrand, V, V_inf)

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


