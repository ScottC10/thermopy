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

        H_r = int[V->inf](T * delP/delT - P)dV + (PV - RT)

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
    V_inf = 100 * V_ideal


    #Integrated over log space to improve accuracy over smaller truncated V_inf
    def log_integrand(x):
        v = np.exp(x)
        return (T * _dpdt(EoS, T, v) - EoS.P_from_TV(T, v))*v

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

            S_r = int[v->inf](delP/delT - R/V)dV

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
    Z = state.Z
    V = state.V

    V_ideal = (R*T) / P
    V_inf = 100 * V_ideal

    #Integrated over log space using substitution u = logv, Jacobian vdV added
    def log_integrand(x):
        v = np.exp(x)
        return (_dpdt(EoS, T, v) - R/v)*v

    integral,_ = quad(
        log_integrand,
        np.log(V_inf),
        np.log(V)
    )

    return integral + R*np.log(Z)



def _dpdt(EoS,T, V):
    h = 1e-5*max(abs(T), 1)

    return (
        EoS.P_from_TV(T+h, V)
         - EoS.P_from_TV(T-h, V)
    ) / (2*h)