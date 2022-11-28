# Test for Injection and envelope
from CoolProp.CoolProp import PropsSI, PhaseSI
from TurboCor import AVB87DA203_50kW, DSG480_4_120kW , DSG480_4R_120kW   # Mitsubishi Kompressor R32
import Fluide

# TODO - Eingangsgrössen für Einspritzung in refstate_1 einbauen. ==> Kompressormodell ändert sich

T_evap = -11.0
T_cond = 55



K0 = 273.15

speed = 100.0  # Kompressordrehzahl in Hz
inj_perc = 0.0  # Einspritzung Massenstrom in %

subcooling = 3.0  # Subcooling nur Kondensator in K
superheat = 5.0  # Überhitzung vor Konpressor (nach IWT)
superheat_inj = 7.0 #  Überhitzung Einspritzung

icomp = 2

if icomp == 1:

    #R32
    refrigerant = Fluide.Refrigerant('R32')    # change it 1
else:
    refrigerant = Fluide.Refrigerant('R1234zeE')     # R1234ze




# 1 Kompressor inlet  [mdot, T, p, h]

p_evap = PropsSI('P', 'T', K0 + T_evap, 'Q', 0.5, refrigerant.name)
h_1 = PropsSI('H', 'T', K0 + T_evap + superheat, 'P', p_evap, refrigerant.name)

refstate_1 = [0, K0 + T_evap + superheat, p_evap, h_1]

p_cond = PropsSI('P', 'T', K0 + T_cond, 'Q', 0.5, refrigerant.name)
h_3 = PropsSI('H', 'T', K0 + T_cond - subcooling, 'P', p_cond, refrigerant.name)

refstate_3 = [0, K0 + T_cond, p_cond, h_3]



# Kompressor

if icomp == 1:
    # Mitsubishi R32 Kompressor
    xc = AVB87DA203_50kW.calculate_fromrpm(speed, T_evap, T_cond)  # xc [ Heizlstg., el Lstg., flow, Strom, Pcool, COP]
else:
    # Danfoss R515B Kompressor
    # xc = DSG480_4_120kW.calculate_fromrpm(58, T_evap, T_cond)  # R515B noch nicht in CoolProp
    xc = DSG480_4R_120kW.calculate_fromrpm(58, T_evap, T_cond)

# refstate_1 Korrigieren - Massenstrom setzen


refstate_1[0] = xc[2] / 3600.0

# Refstate 2

T3E = (T_cond + T_evap) * 0.5 + K0 + superheat_inj
refstate_1E = [inj_perc / 100.0 * xc[2] / 3600.0, T3E, (p_evap + p_cond) * 0.5,
               PropsSI('H', 'T', T3E, 'P', (p_evap + p_cond) * 0.5, refrigerant.name)]

if icomp == 1:
    h2 = AVB87DA203_50kW.h2a(xc[1], refstate_1[3], xc[2] / 3600.0, refstate_1E[0], refstate_1E[3])
else:
    h2 = refstate_1[3] + xc[1] / (xc[2] / 3600.0)   #fuer FixScrolls

thg = PropsSI('T', 'H', h2, 'P', refstate_3[2], refrigerant.name)

if icomp == 1:
    thg_comp = 0.0
else:
    thg_comp = xc[6]

refstate_2 = [xc[2] / 3600.0, thg, refstate_3[2], PropsSI('H', 'T', thg, 'P', refstate_3[2], refrigerant.name)]


# Isentropic compression

h1isen = h_1
s_isen = PropsSI('S', 'T', refstate_1[1], 'P', refstate_1[2], refrigerant.name)

h2isen = PropsSI('H', 'S', s_isen, 'P', refstate_2[2], refrigerant.name)
thgisen = PropsSI('T', 'H', h2isen, 'P', refstate_2[2], refrigerant.name)

print('THG Isentrop: ', thgisen - K0, 'THG Compressor: ', thg_comp, 'THG calculated: ', thg - K0)

P_isentrop = refstate_1[0] * (h2isen - h1isen)
eta_isentrop = P_isentrop / xc[1]



result = [inj_perc, xc[0], xc[1], xc[5], xc[2], thg - K0, eta_isentrop, P_isentrop, refstate_1E[2] / 1.0e6]

print('Perc. Injection: ', inj_perc, 'Ref. Flow: ', xc[2],  'Inj. Pressure: ', refstate_1E[2] / 1.0e6)

print('Heating P: ', xc[0], 'El. Power: ', xc[1], 'COP: ', xc[5], 'Ref. Flow: ', xc[2], 'Eta is: ', eta_isentrop, 'P_isentrop: ', P_isentrop)

