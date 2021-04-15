# summarized actors
import Fluids
import math


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
