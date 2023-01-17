import math

import numpy as np
import pandas as pd
import os
import time
from scipy.interpolate import griddata
from scipy.optimize import fsolve

import Fluids

import HeatExchanger


class LuftKM_OEMWT:
    def __init__(self):
        print("initiated LuftKM_OEMWT")

    def solvebalance(self, iwtact: HeatExchanger.iwt, air: Fluids.HAir, T2phase, power, Tcond, Tevap, subcool):
        thermoparams = (iwtact, air, power, Tcond, Tevap, subcool)
        # print(" Thermoparams: ", thermoparams)

        # print("Thermo Balance: ", self.balancehx(inflow.temperature, *thermoparams))

        x = fsolve(self.balancehx, T2phase, args=thermoparams)  # inflow.temperature Startwert.

        print(x[0])
        self.twophase_temperature = x[0]

        return x

    def balancehx(self, t2phase, *thermoparams):  # ersetzt
        # print("In balancehx: ", thermoparams)

        iwt, air, power, Tcond, Tevap, subcool = thermoparams
        # power = refstate_1[0] * refstate_1[3] - refstate_3[0] * refstate_3[3]
        # iwt = thermoparams[1]
        # evap = thermoparams[0]
        # air = themoparams[2]

        Tair = air.temperature
        Airvolume = air.vdot
        #        power_iwt = iwt.calc(refstate_3E, refstate_1, REF.name)
        power_iwt = iwt.power
        power_airhx = self.get_power(t2phase, Tcond, Tevap, Tair, Airvolume,
                                     subcool)
        print("In balancehx Evapiwt; ", power_iwt, power_airhx, t2phase, power - (power_airhx + power_iwt))
        return power - (power_airhx + power_iwt)


class Sierra_Evap(LuftKM_OEMWT):
    def __init__(self, fname: str):
        # SH5 5° Überhitzung - ohne IHX
        # Verdampfer und Kondensator Daten getrennt - ebenso mit und ohne IHX - wird sonst zu kompliziert und griddata rechnet zu lange
        super().__init__()
        dname = "..\\Datenfiles\\Sierra\\" + fname
        fcsv = dname + '.csv'
        print(os.getcwd())
        pddata = pd.read_csv(fcsv, sep=';')
        linedat = pddata.to_numpy(copy=True)
        self.setData = linedat[3:, 1:5]  # Schnittstelle 4 Kondensationstemp. Verdampfungstemp. Lufttemp. Luftvolumen
        self.tcdata = linedat[3:, 1:2]
        self.tevapdata = linedat[3:, 2:3]
        self.tairdata = linedat[3:, 3:4]

        self.airvolData = linedat[3:, 4:5]  # Volume airflow copy of 24
        self.airvelData = linedat[3:, 24:25]  # Volume airflow

        self.superheatData = linedat[3:, 37:38]  # Volume airflow
        self.subcoolData = linedat[3:, 39:40]  # Volume airflow
        self.CondPowerSensible = linedat[3:, 47:48]  # Leistung trocken - Spalte 48
        self.Condpressuredrop_Air = linedat[3:, 25:26]  # Druckverlust Luft - Spalte 25
        self.CondTout_Air = linedat[3:, 30:31]  # Dry bulb outlet temperature Air

        self.Condpressuredrop_Ref = linedat[3:, 41:42]  # Druckverlust Refrigerant - Spalte 43
        self.Condcharge_Ref = linedat[3:, 42:43]  # Inhalt Refrigerant - Spalte 44
        self.Condflow_Ref = linedat[3:, 40:41]  # Massenstrom Refrigerant - Spalte 42
        self.kA = linedat[3:, 54:55]  # TODO: evtl. berechnen - aus eps = 1 - exp(-kA/mcp) kA = mcp * -LN(1 - q / qmax)

        self.airResData = linedat[3:, 23:33]
        self.refResData = linedat[3:, 36:43]

    def get_power(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume)
        power = float(griddata(self.setData, self.CondPowerSensible, qinp,
                               method='linear'))
        if math.isnan(power):
            kAact = float(griddata(self.setData, self.kA, qinp,
                                   method='nearest'))
            mcp = Airvolume / 3600.0 * 101325 / 287 / (273.15 + Tair) * 1006.0
            pmax = mcp * (Tair - t2phase)
            eps = 1 - math.exp(-kAact / mcp)
            power = eps * pmax
        return -1 * power  # per Definition Verdampfer negativ

    def get_RefDp(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume)
        value = float(griddata(self.setData, self.Condpressuredrop_Ref, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condpressuredrop_Ref, qinp,
                                   method='nearest'))
        return value

    def get_AirDp(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume)
        value = float(griddata(self.setData, self.Condpressuredrop_Air, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condpressuredrop_Air, qinp,
                                   method='nearest'))
        return value

    def get_Refcharge(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume)
        value = float(griddata(self.setData, self.Condcharge_Ref, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condcharge_Ref, qinp,
                                   method='nearest'))
        return value

    def get_AirTout(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume)
        value = float(griddata(self.setData, self.CondTout_Air, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.CondTout_Air, qinp,
                                   method='nearest'))
        return value


class Sierra_Evap_IHX(LuftKM_OEMWT):
    def __init__(self, fname: str):
        # SH0 0° Überhitzung - mit Nachverdampfung im IHX
        # Verdampfer und Kondensator Daten getrennt - ebenso mit und ohne IHX - wird sonst zu kompliziert und griddata rechnet zu lange
        dname = "..\\Datenfiles\\Sierra\\" + fname
        fcsv = dname + '.csv'
        print(os.getcwd())
        pddata = pd.read_csv(fcsv, sep=';')
        linedat = pddata.to_numpy(copy=True)
        self.setData = linedat[3:,
                       1:6]  # Schnittstelle 5 Kondensationstemp. Verdampfungstemp. Lufttemp. Luftvolumen subcooling
        self.tcdata = linedat[3:, 1:2]
        self.tevapdata = linedat[3:, 3:4]
        self.tairdata = linedat[3:, 4:5]

        self.airvolData = linedat[3:, 5:6]  # Volume airflow copy of 24
        self.airvelData = linedat[3:, 25:26]  # Volume airflow

        self.superheatData = linedat[3:, 38:39]  # Volume airflow
        self.subcoolData = linedat[3:, 2:3]  # Volume airflow
        self.CondPowerSensible = linedat[3:, 48:49]  # Leistung trocken - Spalte 48
        self.Condpressuredrop_Air = linedat[3:, 26:27]  # Druckverlust Luft - Spalte 25
        self.CondTout_Air = linedat[3:, 31:32]  # Dry bulb outlet temperature Air

        self.Condpressuredrop_Ref = linedat[3:, 42:43]  # Druckverlust Refrigerant - Spalte 43
        self.Condcharge_Ref = linedat[3:, 43:44]  # Inhalt Refrigerant - Spalte 44
        self.Condflow_Ref = linedat[3:, 41:42]  # Massenstrom Refrigerant - Spalte 42
        self.kA = linedat[3:, 55:56]

        self.airResData = linedat[3:, 24:34]
        self.refResData = linedat[3:, 37:44]

    def get_power(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume, subcool)
        power = float(griddata(self.setData, self.CondPowerSensible, qinp,
                               method='linear'))
        if math.isnan(power):  # Rückfalllösung NTU kA in Spalte 55 (AV
            kAact = float(griddata(self.setData, self.kA, qinp,
                                   method='nearest'))
            mcp = Airvolume / 3600.0 * 101325 / 287 / (273.15 + Tair) * 1006.0
            pmax = mcp * (Tair - t2phase)
            eps = 1 - math.exp(-kAact / mcp)
            power = eps * pmax
            print(" Sierra_Evap_IHX.get_power: ", t2phase, Tcond, Tair, Airvolume, eps, kAact, -power)
        return -1 * power  # per Definition Verdampfer negativ

    def get_RefDp(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume, subcool)
        value = float(griddata(self.setData, self.Condpressuredrop_Ref, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condpressuredrop_Ref, qinp,
                                   method='nearest'))
        return value

    def get_AirDp(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume, subcool)
        value = float(griddata(self.setData, self.Condpressuredrop_Air, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condpressuredrop_Air, qinp,
                                   method='nearest'))
        return value

    def get_Refcharge(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume, subcool)
        value = float(griddata(self.setData, self.Condcharge_Ref, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condcharge_Ref, qinp,
                                   method='nearest'))
        return value

    def get_AirTout(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool):
        qinp = (Tcond, t2phase, Tair, Airvolume, subcool)
        value = float(griddata(self.setData, self.CondTout_Air, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.CondTout_Air, qinp,
                                   method='nearest'))
        return value


class Sierra_Cond(LuftKM_OEMWT):
    def __init__(self, fname: str):
        # Kondensator keine Abhängigkeit von Tevap und daher superheat und subcool - wird in Kreisprozess berücksichtigt        # SH5 5° Überhitzung - ohne IHX
        # Verdampfer und Kondensator Daten getrennt
        super().__init__()
        dname = "..\\Datenfiles\\Sierra\\" + fname
        fcsv = dname + '.csv'
        print(os.getcwd())
        pddata = pd.read_csv(fcsv, sep=';')
        linedat = pddata.to_numpy(copy=True)
        self.setData = linedat[3:, 1:5]  # Schnittstelle 3 Kondensationstemp.  CompressorDischarge Lufttemp. Luftvolumen
        # self.CondPowerSensible = linedat[3:, 43:44]  # Leistung trocken - Spalte 49
        self.tcdata = linedat[3:, 1:2]
        self.tairdata = linedat[3:, 3:4]
        self.tdischarge = linedat[3:, 2:3]
        self.airvolData = linedat[3:, 4:5]  # Volume airflow copy of 24
        self.airvelData = linedat[3:, 24:25]  # Volume airflow

        self.subcoolData = linedat[3:, 38:39]  #
        self.CondPowerSensible = linedat[3:, 44:45]  # Leistung trocken - Spalte 48
        self.Condpressuredrop_Air = linedat[3:, 25:26]  # Druckverlust Luft - Spalte 25
        self.CondTout_Air = linedat[3:, 30:31]  # Dry bulb outlet temperature Air
        self.Condpressuredrop_Ref = linedat[3:, 40:41]  # Druckverlust Refrigerant - Spalte 43
        self.Condcharge_Ref = linedat[3:, 41:42]  # Inhalt Refrigerant - Spalte 44
        self.Condflow_Ref = linedat[3:, 39:40]  # Massenstrom Refrigerant - Spalte 42
        self.kA = linedat[3:, 52:53]  # Berechnetes kA

        self.airResData = linedat[3:, 23:33]
        self.refResData = linedat[3:, 35:42]  # noch zu prüfen

    def solvebalance(self, iwtact: HeatExchanger.iwt, air: Fluids.HAir, T2phase, power, Tcond, Tevap, subcool, tdisch):
        thermoparams = (iwtact, air, power, Tcond, Tevap, subcool, tdisch)
        # print(" Thermoparams: ", thermoparams)

        # print("Thermo Balance: ", self.balancehx(inflow.temperature, *thermoparams))

        x = fsolve(self.balancehx, T2phase, args=thermoparams)  # inflow.temperature Startwert.

        print(x[0])
        self.twophase_temperature = x[0]

        return x

    def balancehx(self, t2phase, *thermoparams):  # ersetzt
        # print("In balancehx: ", thermoparams)

        iwt, air, power, Tcond, Tevap, subcool, tdisch = thermoparams
        # power = refstate_1[0] * refstate_1[3] - refstate_3[0] * refstate_3[3]
        # iwt = thermoparams[1]
        # evap = thermoparams[0]
        # air = themoparams[2]

        Tair = air.temperature
        Airvolume = air.vdot
        #        power_iwt = iwt.calc(refstate_3E, refstate_1, REF.name)
        power_iwt = iwt.power
        power_airhx = self.get_power(t2phase, Tcond, Tevap, Tair, Airvolume,
                                     subcool, tdisch)
        print("In balancehx Condiwt; ", power_iwt, power_airhx, t2phase, power - (power_airhx + power_iwt))
        return power - (power_airhx + power_iwt)

    def get_power(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool, tdisch):
        qinp = (t2phase, tdisch, Tair, Airvolume)
        power = float(griddata(self.setData, self.CondPowerSensible, qinp,
                               method='linear'))
        if math.isnan(power):  # Rückfalllösung NTU kA in Spalte 55 (AV
            kAact = float(griddata(self.setData, self.kA, qinp,
                                   method='nearest'))
            mcp = Airvolume / 3600.0 * 101325 / 287 / (273.15 + Tair) * 1006.0
            pmax = mcp * (t2phase - Tair)
            eps = 1 - math.exp(-kAact / mcp)
            power = eps * pmax
            print(" NTU: ", t2phase, Tcond, Tair, Airvolume, eps, kAact, power)

        if math.isnan(power):
            kAact = float(griddata(self.setData, self.CondPowerSensible, qinp,
                                   method='nearest'))
        return power

    def get_RefDp(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool, tdisch):
        qinp = (t2phase, tdisch, Tair, Airvolume)
        value = float(griddata(self.setData, self.Condpressuredrop_Ref, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condpressuredrop_Ref, qinp,
                                   method='nearest'))
        return value

    def get_AirDp(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool, tdisch):
        qinp = (t2phase, tdisch, Tair, Airvolume)
        value = float(griddata(self.setData, self.Condpressuredrop_Air, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condpressuredrop_Air, qinp,
                                   method='nearest'))
        return value

    def get_Refcharge(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool, tdisch):
        qinp = (t2phase, tdisch, Tair, Airvolume)
        value = float(griddata(self.setData, self.Condcharge_Ref, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.Condcharge_Ref, qinp,
                                   method='nearest'))
        return value

    def get_AirTout(self, t2phase, Tcond, Tevap, Tair, Airvolume, subcool, tdisch):
        qinp = (t2phase, tdisch, Tair, Airvolume)
        value = float(griddata(self.setData, self.CondTout_Air, qinp,
                               method='linear'))
        if math.isnan(value):
            value = float(griddata(self.setData, self.CondTout_Air, qinp,
                                   method='nearest'))
        return value


name1 = "CoilStudio Belaria 65 20L Modellerstellung Verdampfer"
Sierra_2519_BelPro60_SH5 = Sierra_Evap(name1)

name2 = "CoilStudio Belaria 65 20L Modellerstellung VerdampferSH0"
Sierra_2519_BelPro60_IHX = Sierra_Evap_IHX(name2)

name3 = "CoilStudio Belaria 65 20L Modellerstellung"
Sierra_2519_BelPro60_Cond = Sierra_Cond(name3)

print(Sierra_2519_BelPro60_SH5.get_power(-10, 57.5, -16, 2, 24000.0, 2))
print(Sierra_2519_BelPro60_SH5.get_power(-9, 57.5, -16, 2, 24000.0, 2))
print(Sierra_2519_BelPro60_SH5.get_power(-8, 57.5, -16, 2, 24000.0, 2))
print(Sierra_2519_BelPro60_SH5.get_power(-7, 57.5, -16, 2, 24000.0, 2))
print(Sierra_2519_BelPro60_SH5.get_power(-6, 57.5, -16, 2, 24000.0, 2))
print(Sierra_2519_BelPro60_SH5.get_power(-5, 57.5, -16, 2, 24000.0, 2))
print(Sierra_2519_BelPro60_SH5.get_power(-4, 57.5, -16, 2, 24000.0, 2))

# Tevap = -16
# Tcond = 56
# Tair = -7
# Airvolume = 24600.0
# Superheat = 0
# Subcool = 1
#
# qinp = (Tcond, Tevap,  Tair, Airvolume)
# print(qinp)
# power = float(griddata(Sierra_2519_BelPro60_SH5.setData, Sierra_2519_BelPro60_SH5.CondPowerSensible, qinp, method='linear'))
#
# print(power)
#
# power = float(griddata(Sierra_2519_BelPro60_SH5.setData, Sierra_2519_BelPro60_SH5.CondPowerSensible, qinp, method='linear'))
#
# print(power)
#
# x = (griddata(Sierra_2519_BelPro60_SH5.setData, Sierra_2519_BelPro60_SH5.airResData, qinp, method='linear'))
#
# print(x)
#
# y = (griddata(Sierra_2519_BelPro60_SH5.setData, Sierra_2519_BelPro60_SH5.refResData, qinp, method='linear'))
#
# print(y)
#
# y = Sierra_2519_BelPro60_SH5.get_power(Tcond,  Tevap,  Tair, Airvolume)
# print(y)

# Evap mit IHX Test

# qinp = (Tcond, Tevap,  Tair, Airvolume, Subcool)
# print(qinp)
# power = float(griddata(Sierra_2519_BelPro60_IHX.setData, Sierra_2519_BelPro60_IHX.CondPowerSensible, qinp, method='linear'))
#
# print(power)
#
# Subcool = 30
# qinp = (Tcond, Tevap,  Tair, Airvolume, Subcool)
#
# power = float(griddata(Sierra_2519_BelPro60_IHX.setData, Sierra_2519_BelPro60_IHX.CondPowerSensible, qinp, method='linear'))
#
# print(power)
#
# x = (griddata(Sierra_2519_BelPro60_IHX.setData, Sierra_2519_BelPro60_IHX.airResData, qinp, method='linear'))
#
# print(x)
#
# y = (griddata(Sierra_2519_BelPro60_IHX.setData, Sierra_2519_BelPro60_IHX.refResData, qinp, method='linear'))
#
# print(y)

# Kondensator Test


# Tair = 20
# qinp = (Tcond,  Tair, Airvolume)
# print(qinp)
# power = float(griddata(Sierra_2519_BelPro60_Cond.setData, Sierra_2519_BelPro60_Cond.CondPowerSensible, qinp, method='linear'))
#
# print(power)
#
# x = (griddata(Sierra_2519_BelPro60_Cond.setData, Sierra_2519_BelPro60_Cond.airResData, qinp, method='linear'))
#
# print(x)
#
# y = (griddata(Sierra_2519_BelPro60_Cond.setData, Sierra_2519_BelPro60_Cond.refResData, qinp, method='linear'))
#
# print(y)


# Advanced Test   # Achtung: dieser Test dauert ca. 1/2 Stunde
# start_time = time.time()
#
# for tc in range(30, 60, 10):
#     for airvol in range(20000, 40000, 5000):
#         for ta in range(-10, 10, 5):
#
#             tev = ta - 9.0
#             subc = 1
#             sh = 5
#
#             q1 = (float(tc), tev, float(ta), float(airvol))
#             p1 = float(griddata(Sierra_2519_BelPro60_SH5.setData, Sierra_2519_BelPro60_SH5.CondPowerSensible, q1,
#                                    method='linear'))
#
#             q2 = (tc, tev, ta, airvol, subc)
#             p2 = float(griddata(Sierra_2519_BelPro60_IHX.setData, Sierra_2519_BelPro60_IHX.CondPowerSensible, q2,
#                                    method='linear'))
#             q2 = (tc, tev, ta, airvol, subc+20)
#             p2c = float(griddata(Sierra_2519_BelPro60_IHX.setData, Sierra_2519_BelPro60_IHX.CondPowerSensible, q2,
#                                    method='linear'))
#
#             # qc = (tc, ta + 30, airvol)
#             # p3 = float(griddata(Sierra_2519_BelPro60_Cond.setData, Sierra_2519_BelPro60_Cond.CondPowerSensible, qc,
#             #                        method='linear'))
#
#             # print((tc, tev, ta, airvol, subc+20, p1, p2, p2c, ta + 30, p3))
#             print((tc, tev, ta, airvol,  p1, p2, p2c))
#
#
#
#
# print("--- %s seconds ---" % (time.time() - start_time))
