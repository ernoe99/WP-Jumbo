from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI, PhaseSI
import math
from scipy.optimize import fsolve

from TurboCor import AVB87DA203_50kW    # Mitsubishi Kompressor R32
from HeatExchanger import AirEvaporator, iwt, noiwt
import HeatExchanger
import Fluids
import Fluide
import PumpsandFans



# TEst Evaporator


# Ab hier Eingabe

K0 = 273.15
p0 = 101325.0

# INPUTS
# Cycle
T_cond = 55.0  # Initiale Kondensationstemperatur = Ziel Wasser Austrittstemperatur
T_evap = -6.0  # initiale Verdampfungstemperatur - wird in der Berechnung final bestimmt

Tair_in = 2.0  # Sensible Lufteintrittstemperatur
rFair = 83.5  # rel. Feuchte Eintrittsluft

speed = 120.0  # Kompressordrehzahl in Hz
nf_perc = 70.0  # Fan Drehzahl in %
inj_perc = 0.0  # Einspritzung Massenstrom in %

subcooling = 3.0  # Subcooling nur Kondensator in K
superheat = 5.0  # Überhitzung vor Konpressor (nach IWT)
superheat_inj = 7.0 #  Überhitzung Einspritzung

Volumen = 16000.0  # Physikalischer Input Luftvolumen - überschreibt nf_perc
DpressureFan = 67.0  # Druckverlust des Fans

# Heizung
T_RL = 50.0  # vorgegebne RL Temperatur
dT = 5.0  # Ziel delta T für Heizung
P_HP = 65.0 # Ziel Leistung Wärmepumpe in kW
water35 = Fluids.Heatwater(T_RL, 1.0)
m32 = P_HP / water35.cp /dT * 1000.0
print(m32)

# Komponenten & Objekte

kA_evap = 7270.0  # k*A Verdampfer aus trockener Abstimmung /
kA_iwt = 0.0  # 200.0

refrigerant = Fluide.Refrigerant('R32')
air2 = Fluids.HAir(Tair_in, Volumen, rFair / 100.0)
water35 = Fluids.Heatwater(T_RL, m32)  # 50 kW 30-35 Grad vol in l/s

# Kondensator
B427Hx140 = HeatExchanger.Condenser_fitted(140, 0.192, 750.0, (1058, 64.855, 0.0))
print(B427Hx140.uacoeffs)

#Verdampfer und iwt
EV = AirEvaporator(kA_evap)  # Sierra AL7 für FM1
iwt1 = iwt(kA_iwt)
# iwt1 = noiwt

# Fan
fan_power = PumpsandFans.EBM_W3G910_KU25.get_power(volume=Volumen, dpressure=DpressureFan, Temperature=Tair_in)

refstate_3E = [0.171, 292.15, 20*p0, 220000.0]
refstate_1 = [0.171, 265.15, 6*p0, 540000.0]
xe = EV.calculate(air2, T_evap, refstate_3E, refstate_1, refrigerant)
print(xe)
print('**************')

xc = AVB87DA203_50kW.calculate_fromrpm(speed, T_evap, T_cond)
print(xc)



T_evap_old = 100.0
T_cond_old = -100.0
P_HPact = P_HP

i = 0


# Loop Beginn

while (abs(T_evap_old - T_evap) + abs(T_cond_old - T_cond)) > 1.0e-6:

    i += 1
    T_evap_old = T_evap
    T_cond_old = T_cond

    # States 1, 3 setzen

    # 1 Kompressor inlet  [mdot, T, p, h]

    p_evap = PropsSI('P', 'T', K0 + T_evap, 'Q', 0.5, refrigerant.name)
    h_1 = PropsSI('H', 'T', K0 + T_evap + superheat, 'P', p_evap, refrigerant.name)

    refstate_1 = [0, K0 + T_evap + superheat, p_evap, h_1]

    p_cond = PropsSI('P', 'T', K0 + T_cond, 'Q', 0.5, refrigerant.name)
    h_3 = PropsSI('H', 'T', K0 + T_cond - subcooling, 'P', p_cond, refrigerant.name)

    refstate_3 = [0, K0 + T_cond, p_cond, h_3]

    # Kompressor

    xc = AVB87DA203_50kW.calculate_fromrpm(speed, T_evap, T_cond)  # xc [ Heizlstg., el Lstg., flow, Strom, Pcool, COP]

    # refstate_1 Korrigieren - Massenstrom setzen

    refstate_1[0] = xc[2] / 3600.0

    # Refstate 2

    T3E = (T_cond + T_evap) * 0.5 + K0 + superheat_inj
    refstate_1E = [inj_perc / 100.0 * xc[2] / 3600.0, T3E, (p_evap + p_cond) * 0.5, PropsSI('H', 'T', T3E, 'P', (p_evap + p_cond) * 0.5, refrigerant.name)]

    h2 = AVB87DA203_50kW.h2(xc[1], refstate_1[3], xc[2] / 3600.0, refstate_1E[0], refstate_1E[3])

    thg = PropsSI('T', 'H', h2, 'P', refstate_3[2], refrigerant.name)
    print(thg - K0)

    refstate_2 = [xc[2] / 3600.0, thg, refstate_3[2], PropsSI('H', 'T', thg, 'P', refstate_3[2], refrigerant.name)]

    # Kondensator

    T_cond_A = B427Hx140.solvebalance(water35, P_HPact)  # Neue Kondensationstemperatur als array
    T_cond = T_cond_A[0]   # Setzen neue Kondesationstemperatur in Celsius

    # State after Condenser

    p_cond = PropsSI('P', 'T', K0 + T_cond, 'Q', 0.5, refrigerant.name)
    h_3 = PropsSI('H', 'T', K0 + T_cond - subcooling, 'P', p_cond, refrigerant.name)
    print(PhaseSI('T', T_cond + K0 - subcooling, 'P', p_cond, refrigerant.name))

    refstate_3 = [refstate_2[0], T_cond + K0, p_cond, h_3]

    # Einspritzung

    i1mp = (1 - inj_perc / 100.0)
    h3E = (refstate_3[3]*i1mp - inj_perc / 100.0 * (refstate_1E[3] - refstate_3[3])) / i1mp
    T3E = PropsSI('T', 'H', h3E, 'P|liquid', p_cond, refrigerant.name)

    refstate_3E = [refstate_3[0] * i1mp, T3E, p_cond, h3E]
    # Korrektur refstate_1
    refstate_1[0] = refstate_3E[0]

    # iwt
    #
    # # Set refstate 1I to saturation
    #
    # refstate_1I[0] = refstate_3E[0]
    # refstate_1I[1] = T_evap + K0
    # refstate_1I[2] = p_evap
    # refstate_1I[3] = PropsSI('H', 'P', refstate_1E[2], 'Q', 1.0, refrigerant)
    #
    #
    # o_iwt = iwt1
    # while (hnew - refstate)
    refstate_4 = [refstate_3E[0], K0 + T_evap, refstate_1[2], refstate_3E[3]]



    # Verdampfer - ohne IWT


    power = refstate_1[0] * refstate_1[3] - refstate_3E[0] * refstate_3E[3]
    powerc = refstate_2[0] * (refstate_2[3] - refstate_3[3])
    powere = powerc - xc[1]


    tresult = EV.solvebalance(iwt1, air2, T_evap, powere, refstate_1, refstate_3E, refrigerant)

    print(tresult)

    # states update

    T_evap = tresult[0]   # Setzen neue Verdampfungstemperatur in Celsius

    print(" Aktuelle Iteration i: ", i, T_evap, T_cond, refstate_1, refstate_2, refstate_3, refstate_3E, refstate_4, refstate_1E)
    print (" Powers Iteration: ", i, powerc, xc[1], powere, (powerc - P_HP*1000.0))

    P_HPact = powerc / 1000.0   # P_HPact ist in kW
    m32 = P_HPact / water35.cp / dT * 1000.0
    water35 = Fluids.Heatwater(T_RL, m32)  # 50 kW 30-35 Grad vol in l/s


print("ENDE")

evapo = EV.calculate(air2, T_evap, refstate_3E, refstate_1, refrigerant)
condenso = B427Hx140.calculate(water35, T_cond)

heat_power = water35.capacity(water35.temperature) * dT
evap_power = evapo[1]


print("Heat Pump", air2.temperature, water35.temperature, T_evap, T_cond, heat_power, evap_power, xc)

print(refstate_1, refstate_1[1]-K0)
print(refstate_2, refstate_2[1]-K0)
print(refstate_3, refstate_3[1]-K0)
print(refstate_3E, refstate_3E[1]-K0)
print(refstate_1E, refstate_1E[1]-K0)
print(refstate_4, refstate_4[1]-K0)

print("COP with Fan, Inverter: ", heat_power / (xc[1] *1.02 + fan_power), heat_power, fan_power, xc[1]*1.02)
print("Energy Balance: ", heat_power - evap_power - xc[1], heat_power, evap_power, xc[1])




