#single file in folder in case of change to JSON database and generalised class for species
from dataclasses import dataclass
from typing import TypedDict, NotRequired

class CpCoeffs(TypedDict):
    a : float
    b: float
    c: float
    d: float

class entry(TypedDict):
    formula : str
    Tc      : float
    Pc      : float
    omega   : float
    MW      : float
    Cp_coeffs : NotRequired[CpCoeffs]


"""
All species in database must include:
    -name                                   : str
    -formula                                : str
    -critical temperature [K]               : float
    -critical pressure [Pa]                 : float
    -pitzer acentric factor omega           : float
    -molecular weight                       : float
    -heat capacity correlation coefficients : dict (Shomate Equation)
"""
SPECIES_DATA: dict[str, entry] = {
    "methane" : {
        "formula" : "CH4",
        "Tc"    : 190.6,
        "Pc"    : 4599000,
        "omega" : 0.011,
        "MW"    : 16.043,
        "Cp_coeffs":{"a":19.25,"b":5.213e-2,"c":1.197e-5,"d":-1.132e-8}
    },
    "water" : {
        "formula" : "H2O",
        "Tc"    : 647.1,
        "Pc"    : 22055000,
        "omega" : 0.345,
        "MW"    : 18.015
    },
    "hydrogen sulphide" :{
        "formula" : "H2S",
        "Tc"      : 373.5,
        "Pc"      : 8963000,
        "omega"   : 0.094,
        "MW"      : 34.082,
        "Cp_coeffs" : {"a":26.88412, "b":18.67809/(1000), "c":3.434203/(1000**2), "d": -3.378702/(1000**3)} #Valid 298<T<1400
    }
}
@dataclass(frozen=True)
class CpCoefficients:
    a : float
    b : float
    c : float
    d : float
@dataclass(frozen=True)
class Species:

    formula : str
    Tc : float
    Pc : float
    omega : float
    MW : float
    Cp_coeffs : CpCoefficients | None



def get_species(species:str) -> Species:
    try:
        cur_species = SPECIES_DATA[species]
    except KeyError:
        keys=[key for key in SPECIES_DATA]
        raise Exception(f"{species} is not an available species. Available species: {keys}")
    if "Cp_coeffs" in cur_species:
        return Species(
            formula=cur_species.get("formula"),
            Tc=cur_species.get("Tc"),
            Pc=cur_species.get("Pc"),
            omega=cur_species.get("omega"),
            MW=cur_species.get("MW"),

            Cp_coeffs=CpCoefficients(
                a=cur_species["Cp_coeffs"]["a"],
                b=cur_species["Cp_coeffs"]["b"],
                c=cur_species["Cp_coeffs"]["c"],
                d=cur_species["Cp_coeffs"]["d"]
            )

        )
    else:
        return Species(
            formula=cur_species.get("formula"),
            Tc=cur_species.get("Tc"),
            Pc=cur_species.get("Pc"),
            omega=cur_species.get("omega"),
            MW=cur_species.get("MW"),
            Cp_coeffs=None
        )

