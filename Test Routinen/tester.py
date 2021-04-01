import pyromat as pm
import math
from scipy.optimize import fsolve


class Fluid:
    def __init__(self, temp):
        self.temperature = temp
        self.density = 1.0
        self.vdot = 18.0
        self.mdot = self.vdot * self.density
        self.cp = 4000.0

    def enthalpy(self):
        return self.temperature * self.mdot * self.cp

    def getvector(self):
        return [self.mdot, self.vdot, self.temperature, self.density, self.enthalpy(), self.cp]

    def capacity(self):
        return self.mdot * self.cp


class Brine(Fluid):
    def __init__(self, temp):
        super().__init__(self)  # Das scheint es zu brauchen- TODO noch kären warum
        self.temperature = temp
        self.cp = 3600.0


class Sink:  # This is the building
    def __init__(self, VL, volume, flow, heatload, actors : Actor, sensors : Sensor, incoming: Fluid):
        self.temperature_request = VL
        self.volume = volume
        self.q_heat = heatload



