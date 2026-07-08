import numpy as np

from thermopy.eos import PengRobinson

import matplotlib.pyplot as plt
"""
"methane" : {
        "formula" : "CH4",
        "Tc"    : 190.56,
        "Pc"    : 4599000,
"""
PR = PengRobinson("methane")
Tc = 190.56
T=0.7*Tc


"""List Tests"""
# Temperature = 300 #[K]
# pressures = [n for n in range(1000,100000, 1000)]
# volumes = []

# for p in pressures:
#     volume = VdW.solve(Temperature, p)
#     compressibility = VdW.compressibility(P=p, T=Temperature, V= volume["vapour"])
#     print(volume)
#
#
# plt.plot(volumes, pressures)
# plt.ylabel("pressure")
# plt.xlabel("molar volume")
# plt.show()
#
b = PR.OMEGA * (8.314 * Tc) / 459900

V_range = np.linspace(1.5*b, 5e-3, 1000)
P_values = [PR.P_from_TV(T, V) for V in V_range]

plt.plot(V_range, P_values)
plt.xlabel("molar volume")
plt.ylabel("pressure")
plt.show()