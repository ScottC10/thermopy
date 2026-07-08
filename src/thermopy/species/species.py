#single file in folder in case of change to JSON database and generalised class for species


"""
All species in database must include:
    -name                                   : str
    -formula                                : str
    -critical temperature [K]               : float
    -critical pressure [Pa]                 : float
    -pitzer acentric factor omega           : float
    -molecular weight                       : float
    -heat capacity correlation coefficients : dict
"""



SPECIES_DATA = {
    "methane" : {
        "formula" : "CH4",
        "Tc"    : 190.6,
        "Pc"    : 4599000,
        "omega" : 0.012,
        "MW"    : 16.043,
        "Cp_coeffs":{"a":19.25,"b":5.213e-2,"c":1.197e-5,"d":-1.132e-8}
    },
    "water" : {
        "formula" : "H2O",
        "Tc"    : 647.1,
        "Pc"    : 22055000,
        "omega" : 0.345,
        "MW"    : 18.015
    }
}