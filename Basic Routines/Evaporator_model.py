from scipy.interpolate import griddata
from scipy.optimize import fsolve
import numpy as np
import os

# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
import math
import pandas as pd

# Beginning the functions
def tdelta (tout):
    return -7 + min(tout - 2, 0)*-1.0/9.0 + max(tout - 2, 0) * -0.2

def ntu(kA, x):
    mdot, p, tair, tevap = x
    cpmin = (mdot * 1006.0)
    ntu = kA / cpmin
    eps = (1 - math.exp(-ntu))
    pmax = cpmin * (tair - tevap)
    return p - eps * pmax

def ntu_calc(kA, x):
    mdot, p, tair, tevap = x
    cpmin = (mdot * 1006.0)
    ntu = kA / cpmin
    eps = (1 - math.exp(-ntu))
    pmax = cpmin * (tair - tevap)
    return eps * pmax

print(ntu(2000.0,[4.0, 23.e3, 2, -5]))

kA0 = np.array(1000.0)

kA = fsolve(ntu, kA0, [4.0, 23.0e3, 2, -5])


# Eigenschaften Evaporator TODO als Oject
area = 1.5 * 1.76  # Flaeche des Verdampfers

vref = 1.69
SHref = 5
Tkref = 55.0


t1 = -7
v1 = 1.7
SH1 = 5.0
Tk1 = 55.0

print(tdelta(t1))
print(tdelta(2))
print(tdelta(7))
print(tdelta(5))



# x = np.array[[v1/vref, tdelta(t1)/tdelta(2), Tk1/Tkref, SH1/SHref],
#     [1.05*v1/vref, tdelta(t1+14)/tdelta(2), Tk1/Tkref, SH1/SHref]]

# y = [1, 2]
# xi = [1.0, 1.0, 1.0, 1.0]
#
# z = griddata(x, y, xi,  method='linear', fill_value=99)
#
# print(x)
# print(z)

print(os.getcwd())

file = open('..\\Datenfiles\\6_AL100_7,94mm_70R_6RR_L1500_H1750_40K.csv', 'r')
headers = ["Tair", "v", "Tk", "Tevap", "SH", "Psensibel" ]
data = pd.read_csv(file, delimiter = ";",header=0,names=headers)


sm = (data.Tair - data.Tevap)
td = sm.copy()
A_evaporator = 1.5 * 1.74

data["Tdelta"] = (data.Tair - data.Tevap)
data["rho"] = 101325.0 /((data.Tair + 273.15) * 287)
data["mdot"] = data.rho * data.v * A_evaporator
# data["kA"] = data.Psensibel
kA = td.copy()

print(sm.size)

for i in range(sm.size):
    #    td[i] = tdelta(sm[i]) / tdelta(data.Tair[i])
    td[i] = data.Tdelta[i] / tdelta(data.Tair[i])
    kA[i] = fsolve(ntu, kA0, [data.mdot[i], data.Psensibel[i], data.Tair[i], data.Tevap[i]])[0]

data["Td_rel"] = td
data["v_rel"] = data.v/vref
data["Tk_rel"] = data.Tk/Tkref
data["SH_rel"] = data.SH / SHref
data["kA"] = kA

# y = data.v/vref
# z = data.Tk/Tkref
# zz = data.SH / SHref
# Ps = data.Psensibel
#
#
#
# #generate new grid
# X, Y, Z, ZZ = np.mgrid[-20:20:40j, 0:2:10j, 0:2:10j, 0:2:10j]
#
# #interpolate "data.v" on new grid "inter_mesh"
# V = griddata((td, y, z, zz), Ps, (X, Y, Z, ZZ), method='nearest')

# erg = sm.copy()
# kA_erg = sm.copy()
#
# for i in range(sm.size):
#     erg[i] = griddata((td, y, z, zz), Ps, (td[i], y[i], z[i], zz[i]), method='linear')
#     kA_erg[i] = griddata((td, y, z, zz), data.kA, (td[i], y[i], z[i], zz[i]), method='linear')
#     print(f"Wert: ", i, "P sensibel:" , erg[i], "kA: ", kA_erg[i])
#
# v1 = np.ones(25)
# v2 = [i/10.0 for i in range(1, 26)]
# v3 = v1.copy()
# v4 = np.zeros(25) + 1
#
# res = griddata((td, y, z, zz), data.kA, (v1, v2, v3, v4), method='linear')
#
# # #Plot values
# fig = plt.figure()
# ax=fig.gca(projection='3d')
# sc=ax.scatter(X, Y, Z, c=V, cmap=plt.hot())
# plt.colorbar(sc)
# plt.show()

def Fit_Evaporator(tair, vair, Tkond, Tevap, SH, pd_data):

    # Verdampfer-Fitting ausgehend von den Da
    kA_init = pd_data.kA[0]

    # Geschwindigkeitsabhängigkeit
    kv = 0.2938 * vair**3 - 1.3465 * vair**2 + 2.089 * vair - 0.1026

    # Abhängigkeit Kondensationstemperatur - Basis = 55
    ktk =1 + (55 - Tkond) * 41 / kA_init

    # Abhaengigkeit Grad dT
    kdt = 1 + ((tair - Tevap) + tdelta(tair)) * 172.5 / kA_init

    # Abhaengigkeit Superheat
    ksh = 1 + (5 - SH) * 149.4 / kA_init

    return kA_init * kv * ktk * kdt * ksh

kAFE = np.zeros(sm.size)
PSFE = np.zeros(sm.size)
Perc = np.zeros(sm.size)

for i in range(sm.size):
    #    td[i] = tdelta(sm[i]) / tdelta(data.Tair[i])    td[i] = tdelta(sm[i]) / tdelta(data.Tair[i])
    kAFE[i] = Fit_Evaporator(data.Tair[i], data.v[i], data.Tk[i], data.Tevap[i], data.SH[i], data)
    PSFE[i] = ntu_calc(kAFE[i], (data.mdot[i], data.Psensibel[i], data.Tair[i], data.Tevap[i]))
    Perc[i] = (data.Psensibel[i] - PSFE[i])/data.Psensibel[i]*100.0

data["kA_approx"] = kAFE
data["Ps_calc"] = PSFE
data["Ps_perc"] = Perc

print(data)



