from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI
import math

from TurboCor import AVB87DA203_50kW
import HeatExchanger
import Fluids
import Fluide

K0 = 273.15
p0 = 101325.0



# INPUTS
T_cond = 35.0
T_evap = -6.0
nc_perc = 100.0
nf_perc = 70.0
inj_perc = 5.0
subcooling = 3.0
superheat = 5.0



r32 = Fluide.Refrigerant('R32')
air2 = Fluids.HAir(2, 16000.0, 0.825)
water35 = Fluids.Heatwater(30.0, 2.38)  # 50 kW 30-35 Grad vol in l/s



T3 = 35.0
SC = 5.0

p3 = PropsSI('P', 'T', T3 + K0, 'Q', 0.5, 'R32')
h3 = PropsSI('H', 'T', T3 - SC + K0, 'P', p3, 'R32')
h3sat = PropsSI('H', 'T', T3 + K0, 'Q', 0, 'R32')

T2 = -6.0
Superheat = 5.0
T1 = T2 + Superheat
p2 = PropsSI('P', 'T', T2 + K0, 'Q', 0.5, 'R32')
h1 = PropsSI('H', 'T', T2 + Superheat + K0, 'P', p2, 'R32')
h1sat = PropsSI('H', 'T', T1 + K0, 'Q', 1, 'R32')

class AirEvaporator:  #special for air and evaporation
    def __init__(self, kA):

        self.kA = kA

    def calculate(self, vol, Tairin, Rin, Tevap, SH, href_in, refrigerant, href_out=0.0):  # Rin in %, vol in m3/h, T in Celsius

        # SH is superheat - just in case href_out is not set - means 0
        # AIR first
        rhoair_in = 1.0 / HAPropsSI('V', 'T', K0+Tairin, 'P', p0, 'R', Rin/100.0)
        hair_in = HAPropsSI('H', 'T', Tairin + K0, 'P', p0, 'R', Rin / 100.0)
        xair_in = HAPropsSI('W', 'T', Tairin + K0, 'P', p0, 'R', Rin / 100.0)
        Dewp = HAPropsSI('Tdp', 'T', Tairin + K0, 'P', p0, 'R', Rin / 100.0)
        mdotair = rhoair_in * vol / 3600.0
        masscapa = mdotair * cpair
        qmax = masscapa * (Tairin - Tevap)
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
        else:   #Correctur der Wärmeleistung des Verdampfers mit Kondensation und Vereisung
            xair_out = HAPropsSI('W', 'T', tair_out_dry, 'P', p0, 'R', 1.0)
            Tair_out = tair_out_dry - K0
            mdot_water = mdotair * (xair_in - xair_out)
            qcond = mdot_water * 2.46e6   # Verdunstungsenergie bei ca 20° lt. energie-lexikon.info
            if Tair_out < 0.0:
                qice = mdot_water * 334E3  # Schmelzwärme Eis aus studimpu-physik.de
            else:
                qice = 0

            qfull = qact + qcond + qice

        rhoair_out = 1.0 / HAPropsSI('V', 'H', hair_out, 'P', p0, 'W', xair_out)


        TC_in = Tevap + SH
        PC_in = PropsSI('P', 'T', Tevap + K0, 'Q', 0.5, refrigerant)
        if href_out == 0:
           hC_in = PropsSI('H', 'T', Tevap + SH + K0, 'P', p2, refrigerant)  # Derived from CoolProp via SH
        else:
            hC_in = PropsSI('H', 'T', Tevap + SH + K0, 'P', p2, refrigerant)
        h1sat = PropsSI('H', 'T', T1 + K0, 'Q', 1, refrigerant)

        mdot_ref = qfull / (hC_in - href_in)
        Vair_out = mdotair / rhoair_out

        return [mdot_ref*3600, qfull, qact, qcond + qice, Vair_out*3600.0, Tair_out, href_out, tair_out_dry - K0, hC_in]


EV = AirEvaporator(7270.0)  #Sierra AL7 für FM1


vol = 16000.0
Tairin = 2.0
Rin = 83.8

Tevap1 = -6
SH1 = 5.0
href_in = h3
speed = 111.0



xe = EV.calculate(vol, Tairin, Rin, Tevap1, SH1, href_in, 'R32', 0.0)
print(xe)
print ('**************')

xc = AVB87DA203_50kW.calculate_fromrpm(speed, Tevap1, T3)
print(xc)

print([xe[0], xc[2], xe[0] - xc[2]])
print([xe[1], xc[4], xe[1] - xc[4]])

# Get Hotgas Temperature

perc_inj = 5.0
mdot_refC = xe[0] * (1 + perc_inj/100.0)

t_inj = (Tevap1 + T3) / 2
t_inj_SH = 7.0

p_inj = (p2 + p3) / 2.0
print(p_inj)

h1_inj = PropsSI('H', 'T', t_inj + t_inj_SH + K0, 'P', p_inj, 'R32')
h1_inj_sat = PropsSI('H', 'T', T1 + K0, 'Q', 1, 'R32')

print(h1_inj)
print(h1_inj_sat)

h1_ges = (xe[6]*100.0 + perc_inj*h1_inj) / (100.0 + perc_inj)
print(h1_ges)

mdotx = mdot_refC / 3600.0
print(mdotx, xc[1])
print(xc[1] / mdotx)

h2_ges = h1_ges + xc[1] / mdotx
print(h2_ges)
t_hotgas = PropsSI('T', 'H', h2_ges, 'P', p3, 'R32')
print(t_hotgas - K0)

B427Hx140 = HeatExchanger.Condenser_fitted(140, 0.192, 750.0, (1058, 64.855, 0.0))
W35_100 = Fluids.Heatwater(30.0, 2.38)

print(B427Hx140.uacoeffs)


Kout = B427Hx140.solvebalance(W35_100, 50.0)

h4_iwt_in = xe[6]

print(Kout)

class iwt:  #special for Refrigerant / refrigerant
    def __init__(self, kA):

        self.kA = kA

    def calculate(self, mrefL, hrefL,  mrefG, hrefG, refrigerant, pcond, pevap):  # L liquid G gas overheated or twophase

        TrefL = PropsSI('T', 'H', hrefL, 'P', pcond, refrigerant)
        TrefG = PropsSI('T', 'H', hrefG, 'P', pevap, refrigerant)

        cpL = PropsSI('CP0MASS', 'H', hrefL, 'P', pcond, refrigerant)
        cpG = PropsSI('CP0MASS', 'H', hrefG, 'P', pevap, refrigerant)

        h11 = PropsSI('H', 'T', TrefL, 'P', pcond, refrigerant)
        h11p1 = PropsSI('H', 'T', TrefL+1, 'P', pcond, refrigerant)
        cpL = h11p1 - h11
        h22 = PropsSI('H', 'T', TrefG, 'P', pevap, refrigerant)
        h22p1 = PropsSI('H', 'T', TrefG+1, 'P', pevap, refrigerant)
        cpG = h22p1 - h22

        ntu = self.kA / min(mrefL * cpL, mrefG * cpG )
        q_max = min(mrefL * cpL, mrefG * cpG ) * (abs(TrefL - TrefG))

        cr = min(mrefL * cpL, mrefG * cpG )/ max(mrefL * cpL, mrefG * cpG )

        if cr > 0.999:  # Sonderfall
            eps = ntu / (1 + ntu)
        else:
            eps = (1 - math.exp(-ntu * (1 - cr))) / (1 - cr * math.exp(-ntu * (1 - cr)))

        return q_max * eps

# test iWT


hk1 = 240000.0
pk1 = p3
mk1 = 760/3600

hk2 = 520000.0
pk2 = p2
mk2 = mk1

iwt1 = iwt(2000.0)

print(iwt1.calculate(mk1, hk1, mk2, hk2, 'R32', pk1, pk2))



