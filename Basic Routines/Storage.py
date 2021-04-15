from Fluids import Fluid, Heatwater
from PumpsandFans import Pump
import pandas as pd


class Storage:
    def __init__(self, inflow: Fluid, outflow: Fluid, pump: Pump, volume, temperature):
        self.Pump = pump
        self.Volume = volume
        self.Inflow = inflow
        self.Outflow = outflow
        self.Temperature = temperature

    def calculate(self, frac):
        self.Inflow.vdot = self.Pump.flow(frac)
        newt = self.Temperature + inflow.vdot * (self.Inflow.temperature - self.Temperature) / self.Volume
        self.Outflow.vdot = self.Inflow.vdot
        self.Outflow.temperature = self.Temperature
        self.Temperature = newt

        return self.Temperature


class Storage3mix(Storage):
    def __init__(self, inflowhp: Fluid, outflowhp: Fluid, inflowRL: Fluid, outflowVL: Fluid, volume, init_temperature):
        self.Hpinflow = inflowhp
        self.Hpoutflow = outflowhp
        self.BuildinRL = inflowRL
        self.BuildoutVL = outflowVL
        self.Volume = volume
        self.m1 = volume / 3.0
        self.m2 = volume / 3.0
        self.m3 = volume / 3.0
        self.temperature_1 = init_temperature
        self.temperature_2 = init_temperature
        self.temperature_3 = init_temperature
        self.seconds = 0
        self.data = [{"Seconds: ": self.seconds, "Temperature 1": self.temperature_1,
                      "Temperature 2": self.temperature_2, "Temperature 3": self.temperature_3,
                      "Volume HP": self.Hpinflow.vdot,
                      "Temperature HP in": self.Hpinflow.temperature, "Temperature HP out": self.Hpoutflow.temperature,
                      "Volume Building": self.BuildoutVL.vdot,
                      "Temperature RL": self.BuildinRL.temperature,
                      "Temperature VL": self.BuildoutVL.temperature}]
        self.i = 0

    def calculatex(self, hptout, trl, hpvol, bvol, dt):
        self.seconds += dt
        self.i += 1
        tm11 = self.temperature_1 + dt / self.m1 * (
                hpvol * (hptout - self.temperature_1) + (bvol - hpvol) * (self.temperature_2 - self.temperature_1))
        tm31 = self.temperature_3 + dt / self.m3 * (
                bvol * (trl - self.temperature_3) + (bvol - hpvol) * (self.temperature_3 - self.temperature_2))
        tm21 = self.temperature_2 + dt / self.m2 * ((bvol - hpvol) *
                                                    (self.temperature_1 + self.temperature_3 - 2 * self.temperature_2))
        j = self.i

        self.data.append({"Seconds: ": self.seconds, "Temperature 1": tm11,
                          "Temperature 2": tm21, "Temperature 3": tm31,
                          "Volume HP": hpvol,
                          "Temperature HP in": self.temperature_3, "Temperature HP out": hptout,
                          "Volume Building": bvol,
                          "Temperature RL": trl,
                          "Temperature VL": self.temperature_1})
        self.temperature_1 = tm11
        self.temperature_2 = tm21
        self.temperature_3 = tm31


hpin = Heatwater(50, 10.0)
hpout = Heatwater(55, 10.0)

bvl = Heatwater(50, 15.0)
brl = Heatwater(50, 15.0)

store = Storage3mix(hpin, hpout, brl, bvl, 900.0, 50.0)

for i in range(1, 111, 1):
    hppower = 150.0
    buidingP = 150.0

    hpvol = 10.0
    bvol = 12.0

    hptout = store.temperature_3 + hppower / 4.180 / hpvol
    trl = store.temperature_1 - buidingP / 4.180 / bvol

    store.calculatex(hptout, trl, hpvol, bvol, 10.0)

df = pd.DataFrame(store.data)
with pd.ExcelWriter('..\\output\\Storage.xlsx') as writer:
    df.to_excel(writer, sheet_name='Storage')

