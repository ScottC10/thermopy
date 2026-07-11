from pygments.styles import vs
import matplotlib.pyplot as plt
from thermopy.constants import R
from scipy.integrate import quad
from scipy import differentiate
import numpy as np
from thermopy.eos import EoSResult, EoS, CubicEoS
from typing import Union
from thermopy.species import get_species

def residual_enthalpy(state : EoSResult, EoS : Union[EoS, CubicEoS], cubicEoS=False):
    '''
    calculates the residual enthalpy of the species at specific conditions using a departure function
    solved numerically unless cubicEoS == True
    :param state: EoSResult object formed by an equation of state
    :param EoS: The equation of state used to generate state
    :param cubicEoS: bool deciding which solve method is used default = False
    :return:
    '''
    if cubicEoS:
        pass

    T = state.T
    P = state.P
    V = state.V
    Z = state.Z

    V_ideal = (R*T) / P
    V_inf = 10000 * V_ideal

    def dpdt(T, V):
        h = 1e-5*max(abs(T), 1)

        return (
            EoS.P_from_TV(T+h, V)
            - EoS.P_from_TV(T-h, V)
        ) / (2*h)

    def log_integrand(x):
        v = np.exp(x)
        return (T * dpdt(T, v) - EoS.P_from_TV(T, v))*v

    integral,_ = quad(
        log_integrand,
        np.log(V_inf),
        np.log(V)
    )
    H_R = integral
    H_R += (P*V - R*T)
    return H_R

def residual_entropy(state : EoSResult, EoS : Union[EoS, CubicEoS], cubicEoS=False):
    '''
        calculates the residual entropy of the species at specific conditions using a departure function
        solved numerically unless cubicEoS == True
        :param state: EoSResult object formed by an equation of state
        :param EoS: The equation of state used to generate state
        :param cubicEoS: bool deciding which solve method is used default = False
        :return:
        '''
    pass

