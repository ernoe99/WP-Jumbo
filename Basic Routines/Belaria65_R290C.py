from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI, PhaseSI
import numpy as np
import pandas as pd
import time
import openpyxl
import xlwt
import math
from scipy.optimize import fsolve

import PumpsandFans
from TurboCor import VZN175  # Emerson Kompressor
import HeatExchanger
import HXfitting
import LuftKM_OEMWT
import Fluids

import Fluide

# Umrechnungskonstanten für CoolProp

K0 = 273.15
p0 = 101325.0

# Aktuelle Arbeitsversion für R290 Danfoss Kompressor OEM Luftwaermetauscher von Sierra für 4 FP063 von Ziehl Abegg
# PWT von Sierra aktuell noch nicht fix.

# BelariaPro65H             Basis WP mit ca 45 Pa Druckverlust - 4 * PF063 von ZA, Verdampfer 4 rows Sierra, Swep PWT, Danfoss Polynome VZN175
# BelariaPro65C             Wie H nur mit PWT für Kühlung
# BelariaPro65H_IHX         Wie BelariaPro65H mit IHX
# BelariaPro65C_IHX         wie H nur mit PWT für Kühlung

# Daten PWT siehe Excel Auslegungspunkte-Gewichtungen/Komponentenmodelle

# Interface uaf (t1 a0 a1(v) a2(v**2) t2 b0 b1(v) b2(v**2) am Ende multipliziert mit der Fläche

B220Hx128 = HeatExchanger.OnephaseHX_fitted(128, 0.10078125, 2000,
                                            (55.0, 1072.0, 1046.2, 0.0, 5.0, 771.82, 819.78, 0.0))
# Trenntauscher = HeatExchanger.OnephaseHX_fitted(128, 0.10078125, 2000, (55.0, 83.1, 81.1, 0.0, 5.0, 59.83, 63.55, 0.0))
# sole = Fluids.Brine(56.3, 2.906)   # 57 kW mit 6 Grad delta T 3600 kJ/kg
# wasser = Fluids.Heatwater(50, 2.726)   # 57 kW mit 5 Grad delta T
#
# print(Trenntauscher.calculate(sole, wasser, 1))
# print(Trenntauscher.uaf(sole)/Trenntauscher.number_of_plates/Trenntauscher.area_of_plate)
#
#
# sole = Fluids.Brine(55.84, 2.163)
# wasser = Fluids.Heatwater(50, 1.913)
#
# print(Trenntauscher.calculate(sole, wasser, 1))
# print(Trenntauscher.uaf(sole)/Trenntauscher.number_of_plates/Trenntauscher.area_of_plate)
#
# sole = Fluids.Brine(55.46, 1.175)   # 57 kW mit 6 Grad delta T 3600 kJ/kg
# wasser = Fluids.Heatwater(50, 1.2)   # 57 kW mit 5 Grad delta T
#
# print(Trenntauscher.calculate(sole, wasser, 1))
# print(Trenntauscher.uaf(sole)/Trenntauscher.number_of_plates/Trenntauscher.area_of_plate)
#
#
# sole = Fluids.Brine(56.66, 4.158)   # 57 kW mit 6 Grad delta T 3600 kJ/kg
# wasser = Fluids.Heatwater(50, 3.7)   # 57 kW mit 5 Grad delta T
#
# print(Trenntauscher.calculate(sole, wasser, 1))
# print(Trenntauscher.uaf(sole)/Trenntauscher.number_of_plates/Trenntauscher.area_of_plate)


B250ASHx90 = HeatExchanger.Condenser_fitted(90, 0.12667, 2000, (57, 616.19, 134.63, 2.54, 778.68, 350.8))
F250ASHx90 = HeatExchanger.Evaporator_fitted(90, 0.12667, 2000, (57, 616.19, 134.63, 2.54, 778.68, 350.8))
ZA_FP063_4 = PumpsandFans.OEMFan("ZA_FP063", 1.182)

HPwater = Fluids.Brine(12.0, 2.851)
Pp = 57.0

print(F250ASHx90.solvebalance(HPwater, -Pp))

# Test

noiwt = HeatExchanger.noiwt


# The classes

class Belariapro65C:
    def __init__(self, cmp, kd, ev, ihx, tt, fan):
        self.IHX = ihx
        self.KD = kd
        self.EV = ev
        self.CMP = cmp
        self.TT = tt
        self.FAN = fan
        self.REF = Fluide.Refrigerant('R290')
        self.noFans = 4

    def operating_point(self, Tair=42.0, Rair=0.1, Vair=24000.0, Tvl=7.0, Trl=12.0, Power=57.0, TTdT=1.0):

        self.Tair = Tair
        self.Rair = Rair
        self.Vair = Vair  # später über Fan rpm lösen
        self.Tvl = Tvl
        self.Trl = Trl
        self.Power = Power * 1000
        self.Vdot = self.Power / 4184.0 / (
                    Trl - Tvl)  # Achtung Wert aus Fluids Heatwater fix, rl und vl für Kühlen umgekehrt
        self.TTdT = -TTdT  # Temperaturdelta im Trenntauscher Eintritt Achtung negativ beim Kühlen

        #   self.Air = Fluids.HAir(self.Tair, self.Vair, self.Rair)  # besser ausserhalb
        self.Heatwater_out = Fluids.Heatwater(self.Tvl, self.Vdot)
        self.Heatwater_in = Fluids.Heatwater(self.Trl, self.Vdot)

        if (self.TT == 0):
            self.HPwater_out = self.Heatwater_out
            self.HPwater_in = self.Heatwater_in
        else:
            self.HPwater_out = Fluids.Brine(self.Heatwater_out.temperature + self.TTdT, self.Heatwater_out.vdot)
            vdot_cal = self.Heatwater_out.cp / self.HPwater_out.cp * self.Heatwater_out.vdot
            self.HPwater_in = Fluids.Brine(self.Heatwater_in.temperature + self.TTdT, vdot_cal)
            self.HPwater_out.vdot = self.HPwater_in.vdot
            TTpower = self.TT.calculate(self.HPwater_out, self.Heatwater_in, 1)  # 1 für Gegenstrom
            # Anpassung an TTpower - noch falsch
            while abs(abs(TTpower) - self.Power) > 10.0:
                self.HPwater_out.temperature = self.Heatwater_in.temperature - self.Power / -TTpower * \
                                               (self.Heatwater_in.temperature - self.HPwater_out.temperature)
                TTpower = self.TT.calculate(self.HPwater_out, self.Heatwater_in, 1)  # 1 für Gegenstrom

            self.HPwater_in.temperature = self.HPwater_out.heat(self.Power)
            self.TTdT = self.Heatwater_out.temperature - self.HPwater_out.temperature
            # set to inflow condition Heatwater_in+TTdT
            # shiftT = self.HPwater_in.temperature - (self.Heatwater_in.temperature + self.TTdT)
            # self.HPwater_in.temperature = (self.Heatwater_in.temperature + self.TTdT)
            # self.HPwater_out.temperature = self.HPwater_out.temperature + shiftT
            TTpower = self.TT.calculate(self.HPwater_out, self.Heatwater_in, 1)  # 1 für Gegenstrom
            print("Power Control im Trenntauscher: ", TTpower, self.Power)

        self.air_in = Fluids.HAir(Tair, Vair, Rair)

        print("Nach Trenntauscher: ", self.HPwater_out.vdot, self.HPwater_in.temperature, self.HPwater_out.temperature)

        # INPUTS
        # Cycle

        T_cond = 57.0  # Initiale Kondensationstemperatur = Ziel Wasser Austrittstemperatur
        T_evap = 2.0  # initiale Verdampfungstemperatur - wird in der Berechnung final bestimmt

    def loopforPower(self, T_cond, T_evap, subcool, superheat, speed):

        T_evap_old = 100.0  # fuer Iteration
        T_cond_old = -100.0
        i = 0

        while (abs(T_evap_old - T_evap) + abs(T_cond_old - T_cond)) > 1.0e-3 and i < 10:

            i += 1
            T_evap_old = T_evap
            T_cond_old = T_cond

            # States 1, 3 setzen

            # 1 Kompressor inlet  [mdot, T, p, h]

            p_evap = PropsSI('P', 'T', K0 + T_evap, 'Q', 0.5, self.REF.name)
            h_1 = PropsSI('H', 'T', K0 + T_evap + superheat, 'P', p_evap, self.REF.name)
            d_1 = PropsSI('D', 'T', K0 + T_evap + superheat, 'P', p_evap, self.REF.name)

            refstate_1 = [0, K0 + T_evap + superheat, p_evap, h_1, d_1]

            p_cond = PropsSI('P', 'T', K0 + T_cond, 'Q', 0.5, self.REF.name)
            h_3 = PropsSI('H', 'T', K0 + T_cond - subcool, 'P', p_cond, self.REF.name)
            d_3 = PropsSI('D', 'T', K0 + T_cond - subcool, 'P', p_cond, self.REF.name)

            refstate_3 = [0, K0 + T_cond - subcool, p_cond, h_3, d_3]

            # Kompressor

            xc = self.CMP.calculate_direct(speed, T_evap,
                                           T_cond)  # xc [ Heizlstg., el Lstg., flow, Strom, Pcool, COP, Tdischarge]

            Pelectric = xc[1]

            # Refstate 2

            kthg = K0 + xc[6]
            h_2 = PropsSI('H', 'T', kthg, 'P', p_cond, self.REF.name)  # Status h2 direkt
            d_2 = PropsSI('D', 'T', kthg, 'P', p_cond, self.REF.name)

            refstate_2 = [xc[2], kthg, refstate_3[2], h_2, d_2, self.REF.name]

            # Kondensator

            #  T_cond_A = self.KD.solvebalance(self.HPwater_in, self.Power)  # Neue Kondensationstemperatur als array
            # Beim Chiller ist die Leistung auf der Verdampferseite anzugeben. Im Kondensator muss die Kompressorleistung mit abgeführt werden
            # Aufruf hier mit xd[6] Discharge Temperature solvebalance überschrieben für Kondensatorbetrieb
            T_cond_A = self.KD.solvebalance(self.IHX, self.air_in, T_cond, self.Power + Pelectric,
                                            T_cond, T_evap, subcool, xc[6])  # Neue Kondensationstemperatur als array
            T_cond = T_cond_A[0]  # Setzen neue Kondesationstemperatur in Celsius
            # State after Condenser

            p_cond = PropsSI('P', 'T', K0 + T_cond, 'Q', 0.5, self.REF.name)
            h_3 = PropsSI('H', 'T', K0 + T_cond - subcool, 'P', p_cond, self.REF.name)
            d_3 = PropsSI('D', 'T', K0 + T_cond - subcool, 'P', p_cond, self.REF.name)
            # print(PhaseSI('T', T_cond + K0 - subcool, 'P', p_cond, self.REF.name))

            refstate_3 = [refstate_2[0], T_cond + K0 - subcool, p_cond, h_3, d_3]
            refstate_3E = refstate_3  # ohne IHX und Einspritzung

            # ohne iwt

            d_4 = PropsSI('D', 'H', refstate_3E[3], 'P', p_cond, self.REF.name)  # hier H statt T wegen 2-phasen Geb.
            refstate_4 = [refstate_3[0], K0 + T_evap, refstate_1[2], refstate_3[3], d_4]
            # self.IHX.power = 0.0

            # Verdampfer - ohne IWT

            power = refstate_3[0] * refstate_1[3] - refstate_3[0] * refstate_3[3]

            #             tresult = self.EV.solvebalance(self.IHX, self.air_in, T_evap, -power, T_cond, T_evap, subcool)
            tresult = self.EV.solvebalance(self.HPwater_in,
                                           -power / 1000.0)  # Hier muss Power durch 1000 dividiert werden - alte HeatExchanger routine

            print(tresult)

            # states update

            T_evap = tresult[0]  # Setzen neue Verdampfungstemperatur in Celsius

            print(" Aktuelle Iteration i:  ", i, T_evap, T_cond, refstate_1, refstate_2, refstate_3, refstate_4)
            print("   Power:  EV, CMP, Power Inp, KD: ", power, Pelectric, self.Power,
                  refstate_2[0] * (refstate_2[3] - refstate_3[3]))
            print(" Abbruch bei: ", i, abs(T_evap_old - T_evap) + abs(T_cond_old - T_cond))

        self.Tcond = T_cond
        self.Tevap = T_evap
        self.Tdischarge = xc[6]
        self.refstate_1 = refstate_1
        self.refstate_2 = refstate_2
        self.refstate_3 = refstate_3
        self.refstate_3E = refstate_3E
        self.refstate_4 = refstate_4
        self.EV_power = power
        self.CMP_power = Pelectric
        self.KD_power = refstate_2[0] * (refstate_2[3] - refstate_3[3])
        self.IHX_power = refstate_2[0] * (refstate_3[3] - refstate_3E[3])
        self.superheat = superheat
        self.subcool = subcool
        self.REF_massflow = refstate_2[0]
        self.dpair = HP.KD.get_AirDp(T_evap, T_cond, T_evap, self.air_in.temperature, self.air_in.vdot, subcool,
                                     self.Tdischarge)
        self.dpref = HP.KD.get_RefDp(T_evap, T_cond, T_evap, self.air_in.temperature, self.air_in.vdot, subcool,
                                     self.Tdischarge)
        self.charge = HP.KD.get_Refcharge(T_cond, T_cond, T_evap, self.air_in.temperature, self.air_in.vdot, subcool,
                                          self.Tdischarge)
        self.tairout = HP.KD.get_AirTout(T_cond, T_cond, T_evap, self.air_in.temperature, self.air_in.vdot, subcool,
                                         self.Tdischarge)


# Setup the Heatpump

Belariapro65C_SH5 = Belariapro65C(VZN175, LuftKM_OEMWT.Sierra_2519_BelPro60_Cond, F250ASHx90, noiwt, B220Hx128,
                                  ZA_FP063_4)
HP = Belariapro65C_SH5

# Setup operating point
vl = 7.0  # gewünschte Vorlauftemperatur - Achtung Kompressor Limits noch nicht betrachtet - nur Polynome
vair = 24000.0  # Setzen des Volumenstroms - später eventuell fan Power in %
tair = 42.0  # Aussentemperatur Achtung: Es wird nur trocken gerechnet
pp = 50.0  # Zielleistung Verdampfer = PHX - sollte in etwa erreichbar sein. Wird nach der ersten Iteration korrigiert.
sp = 120.0  # Kompressor Drehzahl Vorgabe in rps

tv = vl - 5.0  # Schätzwert für die Verdampfungstemperatur - wird korrigiert hier auf Wasser bezogen

rl = vl + 5.0  # Vorgabe der Spreizung im Heizwasser
tc = tair + 12  # Schätzwert für die initiale Kondensationstemperatur - Delta = 12.0 Grad

# Belariapro65H_SH5.operating_point(Tair=2, Rair=0.1, Vair=24000.0, Tvl=35.0, Trl=30.0, Power=77.0, TTdT=1.0)
# Belariapro65H_SH5.operating_point(Tair=-7, Rair=0.1, Vair=24000.0, Tvl=55.0, Trl=50, Power=55.3, TTdT=1.0)
Belariapro65C_SH5.operating_point(Tair=tair, Rair=0.1, Vair=vair, Tvl=vl, Trl=rl, Power=pp, TTdT=1.0)

# Belariapro65H_SH5.loopforPower(T_cond=37.5, T_evap=-16.0, subcool=2, superheat=2, speed=140)  # Tkond, Tevap, subcool, superheat, speed
# Belariapro65H_SH5.loopforPower(T_cond=57.5, T_evap=-16.0, subcool=2, superheat=5, speed=140)  # Tkond, Tevap, subcool, superheat, speed
Belariapro65C_SH5.loopforPower(T_cond=tc, T_evap=tv, subcool=2, superheat=5,
                               speed=sp)  # Tkond, Tevap, subcool, superheat, speed

print(Belariapro65C_SH5.FAN.get_power(vair, 45.0, tair) * 4)
pf = Belariapro65C_SH5.FAN.get_power(vair, 45.0, tair) * 4

# Ab hier die Loop aus dem Excel - Setpoints.xlsx

# Read data from "..//INPUT//HPSetpoints.xlsx"

setpoints = pd.read_excel("..//INPUT//CHSetpoints.xlsx")

# Loop over setpoints

colnam_hpin = ["HPFin vol", "HPFin T", "HPFout T"]  # 3
colnam_powers = ["Pheat", "PComp", "Pevap", "Pihx", "EER", "PFan", "EER sum"]  # 7
colnam_Ref = ["Tevap", "Tcond", "Massflow Ref", "Evap dp Ref", "Evap Charge Ref", "Superheat", "Subcool"]  # 5
colnam_Air = ["Evap dp Air", "Evap T out"]  # 2

colnams = colnam_hpin + colnam_powers + colnam_Ref + colnam_Air

hpinx = pd.DataFrame([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                     columns=colnams)
Result_df = pd.concat([setpoints, hpinx], axis=1)  # axis=1 entlang der Spalten

headers = ["mflow 1", "Temperature 1", "Pressure 1", "Enhalpy 1", " Density 1",
           "mflow 2", "Temperature 2", "Pressure 2", "Enhalpy 2", " Density 2",
           "mflow 3", "Temperature 3", "Pressure 3", "Enhalpy 3", " Density 3",
           "mflow 3E", "Temperature 3E", "Pressure 3E", "Enhalpy 3E", " Density 3E",
           "mflow 4", "Temperature 4", "Pressure 4", "Enhalpy 4", " Density 4"]

d1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
hpinx1 = pd.DataFrame([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], columns=headers)

Result_df = pd.concat([Result_df, hpinx1], axis=1)

# print(Result_df)
start_time = time.time()

for i in setpoints.index:

    val = setpoints.loc[i]
    rl = val["T Vorlauf"] + val["dT Wasser"]
    Belariapro65C_SH5.operating_point(Tair=val["T Luft ein"], Rair=0.1, Vair=val["Vol Luft [m3/h]"],
                                      Tvl=val["T Vorlauf"],
                                      Trl=rl, Power=val["Kühlleistung [kW]"], TTdT=1.0)


    Belariapro65C_SH5.loopforPower(T_cond=val["T Luft ein"] + 12.0, T_evap=val["T Vorlauf"] - 5.0, subcool=2,
                                   superheat=2,
                                   speed=val[
                                       "Speed Cmp [rps]"])  # Tkond, Tevap, subcool, superheat, speed initialisiert

    # Korrektur, falls für Leistung

    if val["Kühlleistung [kW]"] > 1.01 * HP.EV_power / 1000.0:

        # Angeforderte Leistung ist grösser als berechnete - Reduktion der angeforderten Leistung -
        # sollte nur bei Pmax Anforderung sein. Reduktion der Leistungsanforderung bei gleicher Kompressordrehzahl

        Belariapro65C_SH5.operating_point(Tair=val["T Luft ein"], Rair=0.1, Vair=val["Vol Luft [m3/h]"],
                                          Tvl=val["T Vorlauf"],
                                          Trl=rl, Power=HP.EV_power / 1000.0, TTdT=HP.TTdT)

        Belariapro65C_SH5.loopforPower(T_cond=HP.Tcond, T_evap=HP.Tevap, subcool=2,
                                       superheat=2,
                                       speed=val[
                                           "Speed Cmp [rps]"])  # Tkond, Tevap, subcool, superheat, speed initialisiert

        Result_df.loc[i, "Kühlleistung [kW]"] = HP.EV_power / 1000.0  # Korrektur der Kühlleistung

    elif val["Kühlleistung [kW]"] < 0.99 * HP.EV_power / 1000.0:

        # Angeforderte Leistung ist kleiner als berechnete - Reduktion der angeforderten Kompressordrehzahl linear -
        # sollte nur bei Teillastpunkten sein. Reduktion der Kompressordrehzahl gleiche Leistungsanforderung

        Belariapro65C_SH5.operating_point(Tair=val["T Luft ein"], Rair=0.1, Vair=val["Vol Luft [m3/h]"],
                                          Tvl=val["T Vorlauf"],
                                          Trl=rl, Power=val["Kühlleistung [kW]"], TTdT=HP.TTdT)
        # neue Drehzahl linear
        newspeed = val["Kühlleistung [kW]"] / (HP.EV_power / 1000.0) * val["Speed Cmp [rps]"]
        Belariapro65C_SH5.loopforPower(T_cond=HP.Tcond, T_evap=HP.Tevap, subcool=2,
                                       superheat=2,
                                       speed=newspeed)  # Tkond, Tevap, subcool, superheat, speed initialisiert
        Result_df.loc[i, "Speed Cmp [rps]"] = newspeed  # Korrektur der Kompressordrehzahl

    # TODO:Fan Vair korrigieren mit aktuellem dpair

    Result_df.loc[i, colnams[0]] = HP.HPwater_in.vdot
    Result_df.loc[i, colnams[1]] = HP.HPwater_in.temperature
    Result_df.loc[i, colnams[2]] = HP.HPwater_out.temperature

    Result_df.loc[i, colnams[3]] = HP.KD_power
    Result_df.loc[i, colnams[4]] = HP.CMP_power
    Result_df.loc[i, colnams[5]] = HP.EV_power
    Result_df.loc[i, colnams[6]] = HP.IHX_power
    Result_df.loc[i, colnams[7]] = HP.EV_power / HP.CMP_power

    fanpower = HP.FAN.get_power(HP.air_in.vdot, (HP.dpair * 1.2 + 10.0), HP.air_in.temperature) * HP.noFans
    Result_df.loc[i, colnams[8]] = fanpower
    Result_df.loc[i, colnams[9]] = HP.EV_power / (HP.CMP_power + fanpower)

    Result_df.loc[i, colnams[10]] = HP.Tevap
    Result_df.loc[i, colnams[11]] = HP.Tcond
    Result_df.loc[i, colnams[12]] = HP.REF_massflow
    Result_df.loc[i, colnams[13]] = HP.dpref
    Result_df.loc[i, colnams[14]] = HP.charge
    Result_df.loc[i, colnams[15]] = HP.superheat
    Result_df.loc[i, colnams[16]] = HP.subcool
    Result_df.loc[i, colnams[17]] = HP.dpair
    Result_df.loc[i, colnams[18]] = HP.tairout

    ref = (HP.refstate_1, HP.refstate_2, HP.refstate_3, HP.refstate_3E, HP.refstate_4)
    for j in range(5):  # Writing refstates

        Result_df.loc[i, headers[j * 5]] = ref[j][0]
        Result_df.loc[i, headers[j * 5 + 1]] = ref[j][1] - K0  # Temperatur in Celsius
        Result_df.loc[i, headers[j * 5 + 2]] = ref[j][2] / 1.0E5  # Pressure in bar
        Result_df.loc[i, headers[j * 5 + 3]] = ref[j][3]  # Enthalpy
        Result_df.loc[i, headers[j * 5 + 4]] = ref[j][4]  # Density

Result_df.to_excel('..//output//Chiller_Result.xlsx')

print("--- %s seconds ---" % (time.time() - start_time))

print("ENDE")
