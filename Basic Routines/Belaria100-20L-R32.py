from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI, PhaseSI
import math
from scipy.optimize import fsolve

from TurboCor import AVB87DA203_50kW    # Mitsubishi Kompressor R32
import HeatExchanger
import Fluids

import Fluide



# The classes


class iwt:  # special for Refrigerant / refrigerant
    def __init__(self, kA):

        self.kA = kA

    def calculate(self, mrefL, hrefL, mrefG, hrefG, refrigerant, pcond, pevap):  # L liquid G gas overheated or twophase

        TrefL = PropsSI('T', 'H', hrefL, 'P', pcond, refrigerant)
        TrefG = PropsSI('T', 'H', hrefG, 'P', pevap, refrigerant)

        cpL = PropsSI('CP0MASS', 'H', hrefL, 'P', pcond, refrigerant)
        cpG = PropsSI('CP0MASS', 'H', hrefG, 'P', pevap, refrigerant)

        h11 = PropsSI('H', 'T', TrefL, 'P', pcond, refrigerant)
        h11p1 = PropsSI('H', 'T', TrefL + 1, 'P', pcond, refrigerant)
        cpL = h11p1 - h11
        h22 = PropsSI('H', 'T', TrefG, 'P', pevap, refrigerant)
        h22p1 = PropsSI('H', 'T', TrefG + 1, 'P', pevap, refrigerant)
        cpG = h22p1 - h22

        ntu = self.kA / min(mrefL * cpL, mrefG * cpG)
        q_max = min(mrefL * cpL, mrefG * cpG) * (abs(TrefL - TrefG))

        cr = min(mrefL * cpL, mrefG * cpG) / max(mrefL * cpL, mrefG * cpG)

        if cr > 0.999:  # Sonderfall
            eps = ntu / (1 + ntu)
        else:
            eps = (1 - math.exp(-ntu * (1 - cr))) / (1 - cr * math.exp(-ntu * (1 - cr)))

        return q_max * eps

    def calc(self, refstate_L, refstate_G, refrigerant):

        cpL = PropsSI('C', 'H', refstate_L[3], 'P', refstate_L[2], refrigerant)
        cpG = PropsSI('C', 'H', refstate_G[3], 'P', refstate_G[2], refrigerant)

        ntu = self.kA / min(refstate_L[0] * cpL, refstate_G[0] * cpG)
        q_max = min(refstate_L[0] * cpL, refstate_G[0] * cpG) * (abs(refstate_L[1] - refstate_G[1]))

        cr = min(refstate_L[0] * cpL, refstate_G[0] * cpG) / max(refstate_L[0] * cpL, refstate_G[0] * cpG)

        if cr > 0.999:  # Sonderfall
            eps = ntu / (1 + ntu)
        else:
            eps = (1 - math.exp(-ntu * (1 - cr))) / (1 - cr * math.exp(-ntu * (1 - cr)))

        return q_max * eps


class noiwt(iwt):
    def __init__(self):
        self.kA = 0.0

    def calculate(self, mk1, hk1, mk2, hk2, ref, pk1, pk2):
        return 0.0




noiwt = noiwt()


class AirEvaporator:  # special for air and evaporation
    def __init__(self, kA):

        self.kA = kA

    #    def calculate(self, vol, Tairin, Rin, Tevap, SH, href_in, refrigerant, href_out=0.0):  # Rin in %, vol in m3/h, T in Celsius

    def calculate(self, air: Fluids.HAir, Tevap, refstate_I: list, refstate_O : list,  refrigerant):
        # SH is superheat - just in case href_out is not set - means 0
        # AIR first
        Tairin = air.Ktemperature
        TevapK = 273.15 + Tevap
        p0 = air.pressure
        Rin = air.rF
        vol = air.vdot

        rhoair_in = 1.0 / HAPropsSI('V', 'T', Tairin, 'P', p0, 'R', Rin)
        hair_in = HAPropsSI('H', 'T', Tairin, 'P', p0, 'R', Rin)
        xair_in = HAPropsSI('W', 'T', Tairin, 'P', p0, 'R', Rin)
        Dewp = HAPropsSI('Tdp', 'T', Tairin, 'P', p0, 'R', Rin)
        cpair = HAPropsSI('C', 'T', Tairin, 'P', p0, 'R', Rin)

        mdotair = rhoair_in * vol / 3600.0
        masscapa = mdotair * cpair
        qmax = masscapa * (Tairin - TevapK)
        NTU = self.kA / masscapa

        eps = 1 - math.exp(-NTU)

        qact = qmax * eps

        hair_out = hair_in - qact / mdotair
        tair_out_dry = HAPropsSI('T', 'H', hair_out, 'P', p0, 'W', xair_in)
        hmin_out = HAPropsSI('H', 'T', Dewp, 'P', p0, 'R', 1.0)

        if hair_out > hmin_out:
            xair_out = xair_in
            Tair_out = HAPropsSI('T', 'H', hair_out, 'P', p0, 'W', xair_in) - K0
            qcond = 0
            qice = 0
        else:  # Correctur der Wärmeleistung des Verdampfers mit Kondensation und Vereisung
            xair_out = HAPropsSI('W', 'T', tair_out_dry, 'P', p0, 'R', 1.0)
            Tair_out = tair_out_dry - K0  # TODO Korrektur mit Isenthalpe ???
            mdot_water = mdotair * (xair_in - xair_out)
            qcond = mdot_water * 2.46e6  # Verdunstungsenergie bei ca 20° lt. energie-lexikon.info
            if Tair_out < 0.0:
                qice = mdot_water * 334E3  # Schmelzwärme Eis aus studimpu-physik.de
            else:
                qice = 0

            qfull = qact + qcond + qice

        rhoair_out = 1.0 / HAPropsSI('V', 'H', hair_out, 'P', p0, 'W', xair_out)

        TC_in = refstate_O[1]
        PC_in = refstate_O[2]
        hC_in = refstate_O[3]
        h1sat = PropsSI('H', 'T', T_evap + K0, 'Q', 1, refrigerant)

        href_in = refstate_I[3]

        mdot_ref = qfull / (hC_in - href_in)
        Vair_out = mdotair / rhoair_out

        return [mdot_ref * 3600, qfull, qact, qcond + qice, Vair_out * 3600.0, Tair_out, tair_out_dry - K0,
                hC_in]

    def calc(self, Tevap, air: Fluids.HAir):
        # nur für hxbalance - liefert Tevap neu
        # AIR first
        Tairin = air.Ktemperature
        p0 = air.pressure
        Rin = air.rF
        vol = air.vdot
        TevapK = K0 + Tevap

        rhoair_in = 1.0 / HAPropsSI('V', 'T', Tairin, 'P', p0, 'R', Rin)
        hair_in = HAPropsSI('H', 'T', Tairin, 'P', p0, 'R', Rin)
        xair_in = HAPropsSI('W', 'T', Tairin, 'P', p0, 'R', Rin)
        Dewp = HAPropsSI('Tdp', 'T', Tairin, 'P', p0, 'R', Rin)
        cpair = HAPropsSI('C', 'T', Tairin, 'P', p0, 'R', Rin)

        mdotair = rhoair_in * vol / 3600.0
        masscapa = mdotair * cpair
        qmax = masscapa * (Tairin - TevapK)
        NTU = self.kA / masscapa

        eps = 1 - math.exp(-NTU)

        qact = qmax * eps

        hair_out = hair_in - qact / mdotair
        tair_out_dry = HAPropsSI('T', 'H', hair_out, 'P', p0, 'W', xair_in)
        hmin_out = HAPropsSI('H', 'T', Dewp, 'P', p0, 'R', 1.0)

        if hair_out > hmin_out:
            xair_out = xair_in
            Tair_out = HAPropsSI('T', 'H', hair_out, 'P', p0, 'W', xair_in) - K0
            qcond = 0
            qice = 0
        else:  # Correctur der Wärmeleistung des Verdampfers mit Kondensation und Vereisung
            xair_out = HAPropsSI('W', 'T', tair_out_dry, 'P', p0, 'R', 1.0)
            Tair_out = tair_out_dry - K0
            mdot_water = mdotair * (xair_in - xair_out)
            qcond = mdot_water * 2.46e6  # Verdunstungsenergie bei ca 20° lt. energie-lexikon.info
            if Tair_out < 0.0:
                qice = mdot_water * 334E3  # Schmelzwärme Eis aus studimpu-physik.de
            else:
                qice = 0

            qfull = qact + qcond + qice

        return qfull

    def solvebalance(self, iwtact : iwt, air : Fluids.HAir,  T2phase, power):
        thermoparams = (iwtact, air, power)
        # print(" Thermoparams: ", thermoparams)

        # print("Thermo Balance: ", self.balancehx(inflow.temperature, *thermoparams))

        x = fsolve(self.balancehx, T2phase, args=thermoparams)  # inflow.temperature Startwert.
        # x = fsolve(KD27.balancehx, 34.0,
        #            args=(self.ua, inflow.capacity(), inflow.temperature, compressor.power_condenser()))

        print(x[0])
        self.twophase_temperature = x[0]

        return x

    def balancehx(self, t2phase, *thermoparams):  # ersetzt
        # print("In balancehx: ", thermoparams)

        iwt, air, power = thermoparams
        # power = refstate_1[0] * refstate_1[3] - refstate_3[0] * refstate_3[3]
        # iwt = thermoparams[1]
        # evap = thermoparams[0]
        # air = themoparams[2]

        power_iwt = iwt.calc(refstate_3E, refstate_1, r32.name)
        power_evap = self.calc(t2phase, air)
        print("In balancehx Evapiwt; ", power_iwt, power_evap, t2phase, power - (power_evap + power_iwt))
        return power - (power_evap + power_iwt)

# TEst Evaporator


# Ab hier Eingabe

K0 = 273.15
p0 = 101325.0

# INPUTS
# Cycle
T_cond = 35.0  # Initiale Kondensationstemperatur = Ziel Wasser Austrittstemperatur
T_evap = -6.0  # initiale Verdampfungstemperatur - wird in der Berechnung final bestimmt

Tair_in = 2.0  # Sensible Lufteintrittstemperatur
rFair = 83.5  # rel. Feuchte Eintrittsluft

speed = 100.0  # Kompressordrehzahl in Hz
nf_perc = 70.0  # Fan Drehzahl in %
inj_perc = 5.0  # Einspritzung Massenstrom in %

subcooling = 3.0  # Subcooling nur Kondensator in K
superheat = 5.0  # Überhitzung vor Konpressor (nach IWT)
superheat_inj = 7.0 #  Überhitzung Einspritzung

Volumen = 16000.0  # Physikalischer Input Luftvolumen - überschreibt nf_perc

# Heizung
T_RL = 30.0  # vorgegebne RL Temperatur
dT = 5.0  # Ziel delta T für Heizung
P_HP = 50.0  # Ziel Leistung Wärmepumpe in kW
m32 = P_HP / 4200.0 / dT * 1000.0
print(m32)

# Komponenten & Objekte

kA_evap = 7270.0  # k*A Verdampfer aus trockener Abstimmung / FM1

r32 = Fluide.Refrigerant('R32')
air2 = Fluids.HAir(Tair_in, Volumen, rFair / 100.0)
water35 = Fluids.Heatwater(T_RL, m32)  # 50 kW 30-35 Grad vol in l/s

# Kondensator
B427Hx140 = HeatExchanger.Condenser_fitted(140, 0.192, 750.0, (1058, 64.855, 0.0))
print(B427Hx140.uacoeffs)

#Verdampfer
EV = AirEvaporator(kA_evap)  # Sierra AL7 für FM1


refstate_3E = [0.171, 292.15, 20*p0, 220000.0]
refstate_1 = [0.171, 265.15, 6*p0, 540000.0]
xe = EV.calculate(air2, T_evap, refstate_3E, refstate_1, r32.name)
print(xe)
print('**************')

xc = AVB87DA203_50kW.calculate_fromrpm(speed, T_evap, T_cond)
print(xc)


# test iWT


hk1 = 240000.0
pk1 = 20*p0
mk1 = 760 / 3600

hk2 = 520000.0
pk2 = 6*p0
mk2 = mk1

iwt1 = iwt(200.0)  # 200 - 2000
# iwt1 = noiwt

print(iwt1.calculate(mk1, hk1, mk2, hk2, 'R32', pk1, pk2))

T_evap_old = 100.0
T_cond_old = -100.0

i = 0


# Loop Beginn

while (abs(T_evap_old - T_evap) + abs(T_cond_old - T_cond)) > 1.0e-3:

    i += 1
    T_evap_old = T_evap
    T_cond_old = T_cond

    # States 1, 3 setzen

    # 1 Kompressor inlet  [mdot, T, p, h]

    p_evap = PropsSI('P', 'T', K0 + T_evap, 'Q', 0.5, r32.name)
    h_1 = PropsSI('H', 'T', K0 + T_evap + superheat, 'P', p_evap, r32.name)

    refstate_1 = [0, K0 + T_evap + superheat, p_evap, h_1]

    p_cond = PropsSI('P', 'T', K0 + T_cond, 'Q', 0.5, r32.name)
    h_3 = PropsSI('H', 'T', K0 + T_cond - subcooling, 'P', p_cond, r32.name)

    refstate_3 = [0, K0 + T_cond, p_cond, h_3]

    # Kompressor

    xc = AVB87DA203_50kW.calculate_fromrpm(speed, T_evap, T_cond)  # xc [ Heizlstg., el Lstg., flow, Strom, Pcool, COP]

    Pelectric = xc[0]

    # Refstate 2

    T3E = (T_cond + T_evap) * 0.5 + K0 + superheat_inj
    refstate_1E = [inj_perc/100.0*xc[2]/3600.0, T3E, (p_evap + p_cond)* 0.5, PropsSI('H', 'T', T3E, 'P', (p_evap + p_cond) * 0.5, r32.name)]

    h2 = AVB87DA203_50kW.h2(xc[1], refstate_1[3], xc[2] / 3600.0, refstate_1E[0], refstate_1E[3])

    thg = PropsSI('T', 'H', h2, 'P', refstate_3[2], r32.name)
    print(thg - K0)

    refstate_2 = [xc[2]/3600.0, thg, refstate_3[2], PropsSI('H', 'T', thg + K0, 'P', refstate_3[2], r32.name)]

    # Kondensator

    T_cond_A = B427Hx140.solvebalance(water35, 50.0)  # Neue Kondensationstemperatur als array
    T_cond = T_cond_A[0]   # Setzen neue Kondesationstemperatur in Celsius

    # State after Condenser

    p_cond = PropsSI('P', 'T', K0 + T_cond, 'Q', 0.5, r32.name)
    h_3 = PropsSI('H', 'T', K0 + T_cond - subcooling, 'P', p_cond, r32.name)
    print(PhaseSI('T', T_cond + K0 - subcooling, 'P', p_cond, r32.name))

    refstate_3 = [refstate_2[0], T_cond + K0, p_cond, h_3]

    # Einspritzung

    i1mp = (1 - inj_perc / 100.0)
    h3E = (refstate_3[3]*i1mp - inj_perc / 100.0 * (refstate_1E[3] - refstate_3[3])) / i1mp
    T3E = PropsSI('T', 'H', h3E, 'P|liquid', p_cond, r32.name)

    refstate_3E = [refstate_3[0] * i1mp, T3E, p_cond, h3E]
    # Korrektur refstate_1
    refstate_1[0] = refstate_3E[0]

    #ohne iwt

    refstate_4 = [refstate_3E[0], K0 + T_evap, refstate_1[2], refstate_3E[3]]



    # Verdampfer - ohne IWT


    power = refstate_1[0] * refstate_1[3] - refstate_3E[0] * refstate_3E[3]

    tresult = EV.solvebalance(iwt1, air2, T_evap, power)

    print(tresult)

    # states update

    T_evap = tresult[0]   # Setzen neue Verdampfungstemperatur in Celsius

    print(" Aktuelle Iteration i: ", i, T_evap, T_cond, refstate_1, refstate_2, refstate_3, refstate_3E, refstate_4, refstate_1E)


print("ENDE")

