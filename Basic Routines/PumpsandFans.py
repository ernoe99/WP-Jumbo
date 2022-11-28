# summarized actors
import math

import numpy as np
import pandas as pd
import os
from scipy.interpolate import griddata
import Fluids


class Pump:
    def __init__(self, maxflow, maxpressure, eff):
        self.max_flow = maxflow
        self.max_pressure = maxpressure  #  in bar
        self.efficiency = eff
        self.electric_power = 0.0

    def flow(self, frac):
        return self.max_flow * frac

    def d_pressure(self, fluid: Fluids.Fluid):
        return fluid.vdot / self.max_flow * self.max_pressure

    def calculate(self, fluid : Fluids.Fluid):
        self.electric_power = self.d_pressure(fluid) * fluid.vdot / self.efficiency
        return self.electric_power


class NoPump(Pump):
    def __init__(self):
        super().__init__(0, 0, 1)

    def d_pressure(self, fluid):
        return 0


class Fan:
    def __init__(self, efficiency, max_volume, dpressure):
        self.efficiency = efficiency
        self.max_volume = max_volume
        self.max_dpressure = dpressure
        self.electric_power = dpressure * max_volume / efficiency

    def d_pressure(self, air):
        return self.max_dpressure * (max(air.vdot / self.max_volume, 1.0))**2

    def calculate(self, air: Fluids.WetAir):
        self.electric_power = self.d_pressure(air) * air.vdot / self.efficiency
        return self.electric_power

    def get_air(self, air: Fluids.WetAir, percent):
        air.vdot = math.sqrt(percent * self.max_volume**2)
        self.calculate(air)
        return air


class FlexFan(Fan):
    def __init__(self, efficiency, max_volume, dpressure):
        super().__init__(efficiency, max_volume, dpressure)
        self.modulated = 1


class OEMFan(Fan):
    def __init__(self, fname : str, rho_base):
        self.rho_base = rho_base    # Basis Dichte der Messungen/Daten
        dname = "..\\Datenfiles\\" + fname  # Vielleicht hier ohne "..\\"
        fcsv = dname + '.csv'
        print(os.getcwd())
        pddata = pd.read_csv(fcsv, sep=';')
        linedat = pddata.to_numpy(copy=True)
        self.fanData = linedat[:, 1:3]
        self.fanPower = linedat[:, 3:4]
        self.fanAcoustic = linedat[:, 5:6]
        self.fanRpm = linedat[:, 4:5]

    def get_power(self, volume, dpressure, Temperature):
        rho = 101325 / 287 / (273.15 + Temperature)
        power = float(griddata(self.fanData, self.fanPower, (volume, dpressure), method='linear'))
        if math.isnan(power):
            power = float(griddata(self.fanData, self.fanPower, (volume, dpressure), method='nearest'))

        return power * rho / self.rho_base

    def get_acoustic(self, volume, dpressure, Temperature):

        power = float(griddata(self.fanData, self.fanAcoustic, (volume, dpressure), method='linear'))
        if math.isnan(power):
            power = float(griddata(self.fanData, self.fanAcoustic, (volume, dpressure), method='nearest'))

        return power

    def get_rpm(self, volume, dpressure, Temperature):

        power = float(griddata(self.fanData, self.fanRpm, (volume, dpressure), method='linear'))
        if math.isnan(power):
            power = float(griddata(self.fanData, self.fanRpm, (volume, dpressure), method='nearest'))

        return power

    def get_all(self, volume, dpressure, Temperature):
        return np.array((self.get_rpm(volume, dpressure, Temperature),
                         self.get_power(volume, dpressure, Temperature),
                         self.get_acoustic(volume, dpressure, Temperature)))

EBM_W3G910_KU25 = OEMFan("EBM-W3G910-KU25", 1.15)

print(EBM_W3G910_KU25.get_power(6000.0, 67.4, 25.0))
print(EBM_W3G910_KU25.get_power(16000.0, 67.4, 0.0))
print(EBM_W3G910_KU25.get_power(16000.0, 67.4, -40.0))

print(EBM_W3G910_KU25.get_acoustic(6000.0, 67.4, 25.0))
print(EBM_W3G910_KU25.get_acoustic(16000.0, 67.4, 0.0))
print(EBM_W3G910_KU25.get_acoustic(16000.0, 67.4, -40.0))

print(EBM_W3G910_KU25.get_rpm(6000.0, 67.4, 25.0))
print(EBM_W3G910_KU25.get_rpm(16000.0, 67.4, 0.0))
print(EBM_W3G910_KU25.get_rpm(16000.0, 67.4, -40.0))

ZA_FP063 = OEMFan("ZA_FP063", 1.182)

print(ZA_FP063.get_rpm(6000, 34.5, 25))
print(ZA_FP063.get_acoustic(6000, 45, -20))
print(ZA_FP063.get_power(6500, 45.0, -20))
print(ZA_FP063.get_power(6500, 45.0, 60))
print(ZA_FP063.get_all(6500, 45.0, 60))



ZA_ZN100 = OEMFan("ZA_ZN100-ZILGG", 1.162)

print(ZA_ZN100.get_rpm(16000, 134.5, 25))
print(ZA_ZN100.get_acoustic(16000, 145, -20))
print(ZA_ZN100.get_power(16500, 145.0, -20))
print(ZA_ZN100.get_power(16500, 145.0, 60))
print(ZA_ZN100.get_all(16500, 115.0, 60))

# class EBM_W3G910-KU25(OEMFan):
#     def __init__(self, rho_base = 1.15):
#
#         self.FanData = np.array(((33630, 0), (29845, 80), (25915, 150), (20280, 220), (29025, 0), (25820, 61), (22440, 112),
#                             (17560, 165), (23900, 0), (21260, 41), (184804, 76), (14460, 112), (18780, 0), (16075, 26),
#                             (14520, 47), (11360, 69)))
#         self.FanPower = np.array((1563,
#                              1938,
#                              2245,
#                              2550,
#                              1004,
#                              1255,
#                              1457,
#                              1663,
#                              561,
#                              701,
#                              814,
#                              929,
#                              272,
#                              340,
#                              395,
#                              451))
#
#     def get_power(self, volume, dpressure, Temperature):
#         rho = 101323 / 287 / (273.15 + Temperature)
#         power = float(griddata(self.FanData, self.FanPower, (volume, dpressure), method='linear'))
#         if math.nan(power):
#             power = griddata(self.FanData, self.FanPower, (volume, dpressure), method='nearest')
#
#         return power * rho / self.rho_baseF


