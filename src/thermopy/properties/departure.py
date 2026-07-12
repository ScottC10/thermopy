from pygments.styles import vs
import matplotlib.pyplot as plt
from thermopy.constants import R
from scipy.integrate import quad
from scipy import differentiate
import numpy as np
from thermopy.eos import EoSResult, EoS, CubicEoS, VanDerWaals
from typing import Union
from thermopy.species import get_species

def residual_enthalpy(state : EoSResult, EoS : Union[EoS, CubicEoS], cubicEoS=False) ->float:
    '''
    calculates the residual enthalpy of the species at specific conditions using a departure function
    solved numerically unless cubicEoS == True

        H_r = int[V->inf](T * delP/delT - P)dV + (PV - RT)

    :param state: EoSResult object formed by an equation of state
    :param EoS: The equation of state used to generate state
    :param cubicEoS: bool deciding which solve method is used default = False
    :return:
    '''


    T = state.T
    P = state.P
    V = state.V
    Z = state.Z

    if cubicEoS:
        a, b = EoS._cubic_parameters(T)
        dadt = EoS._da_dT(T)

        if isinstance(EoS, VanDerWaals):
            term1 = -a/V
        else:
            term1 = ( (a - T*dadt) / (b*(EoS.sigma - EoS.epsilon)) ) * np.log((V + EoS.epsilon*b)/(V+ EoS.sigma*b))
        return term1 + (P*V - R*T)

    V_ideal = (R*T) / P
    V_inf = 500 * V_ideal

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

def residual_entropy(state : EoSResult, EoS : Union[EoS, CubicEoS], cubicEoS=False) -> float:
    '''
        calculates the residual entropy of the species at specific conditions using a departure function

            S_r = int[v->inf](delP/delT - R/V)dV

        solved numerically unless cubicEoS == True
        :param state: EoSResult object formed by an equation of state
        :param EoS: The equation of state used to generate state
        :param cubicEoS: bool deciding which solve method is used default = False
        :return:
        '''


    T = state.T
    P = state.P
    Z = state.Z
    V = state.V

    if cubicEoS:


        a, b = EoS._cubic_parameters(T)
        dadt = EoS._da_dT(T)
        if isinstance(EoS, VanDerWaals):
            return -R*np.log(V/(V-b)) + R*np.log(Z)

        term1 = -R*np.log(V/(V-b))
        term2 = -dadt/(b*(EoS.sigma-EoS.epsilon)) *np.log((V+EoS.epsilon*b)/(V+EoS.sigma*b))

        return term1 + term2 + R*np.log(Z)




    V_ideal = (R*T) / P
    V_inf = 500 * V_ideal

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

def residual_gibbs(state : EoSResult, EoS : Union[EoS, CubicEoS], cubicEoS = False) -> float:
    '''

    Uses residual gibbs definition G^R = H^R - TS^R, if EoS used is cubic, analytical result used.

    See .fugacity for equivalent residual gibbs from fugacity coefficient


    :param state: EoSResult object
    :param EoS: EoS used to generate EoSResult
    :return: The residual gibbs free energy
    '''
    T = state.T
    is_cubic = cubicEoS
    h = residual_enthalpy(state, EoS, cubicEoS=is_cubic)
    s = residual_entropy(state, EoS, cubicEoS=is_cubic)

    return h - T*s




def _dpdt(EoS,T, V):
    h = 1e-5*max(abs(T), 1)

    return (
        EoS.P_from_TV(T+h, V)
         - EoS.P_from_TV(T-h, V)
    ) / (2*h)