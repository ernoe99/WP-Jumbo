import math
import pandas as pd

from scipy.optimize import fsolve

from Fluids import Fluid, Heatwater, Brine, Air, WetAir, Burned_Air
from PumpsandFans import Fan


class HeatExchanger:
    def __init__(self, no_plates, area_plate, htc):
        self.number_of_plates = no_plates
        self.area_of_plate = area_plate
        self.heat_transfer_coeff = htc
        self.ua = self.heat_transfer_coeff * self.number_of_plates * self.area_of_plate

    def ntu(self, fluid):
        """

        :type fluid: Fluid
        """
        return self.ua / fluid.capacity(fluid.temperature)

    def calculate(self, hotflow: Fluid, coldflow: Fluid, iflow):
        # ntu = max(self.ntu(hotflow), self.ntu(coldflow))
        print("NTU: ", self.ntu(hotflow), self.ntu(coldflow), self.ua, hotflow.capacity(hotflow.temperature),
              coldflow.capacity(coldflow.temperature))
        ntu = self.ua / min(hotflow.capacity(hotflow.temperature), coldflow.capacity(coldflow.temperature))
        q_max = min(hotflow.capacity(hotflow.temperature), coldflow.capacity(coldflow.temperature)) * \
                (hotflow.temperature - coldflow.temperature)  # max. Transfer in W
        cr = min(hotflow.capacity(hotflow.temperature), coldflow.capacity(coldflow.temperature)) / \
             max(hotflow.capacity(hotflow.temperature), coldflow.capacity(coldflow.temperature))
        print("NTU: ", ntu, cr, q_max, hotflow.temperature, coldflow.temperature)
        eps = 0

        if iflow == 1:  # Gegenstrom
            if cr > 0.999:  # Sonderfall
                eps = ntu / (1 + ntu)
            else:
                eps = (1 - math.exp(-ntu * (1 - cr))) / (1 - cr * math.exp(-ntu * (1 - cr)))
        if iflow == 2:  # Kreuzstrom
            #  eps = 1 - math.exp(qmax**0.22/cr * (math.exp(-cr * qmax**0.78) - 1))
            eps = 1 - math.exp((math.exp(-cr * ntu ** 0.78) - 1) / (cr * ntu ** -0.22))

            print((cr * ntu ** -0.22))
            print(math.exp(-cr * ntu ** 0.78) - 1)
            print(math.exp(-cr * ntu ** 0.78))
            print(-cr)
            print(ntu ** 0.78)
            print(math.exp((math.exp(-cr * ntu ** 0.78) - 1) / (cr * ntu ** -0.22)))

        print("EPS: ", eps, q_max * eps)
        return q_max * eps



class TwophaseHX(HeatExchanger):

    def uaf(self, vol):
        return (self.uacoeffs[0] + self.uacoeffs[1] * vol + self.uacoeffs[2] * vol * vol) * self.area_of_plate * self.number_of_plates

    def ntu(self, fluid: Fluid):
        return self.uaf(fluid.vdot) / fluid.capacity(fluid.temperature)

    def calculate(self, fluid: Fluid, t2phase):
        eps = 1 - math.exp(-self.ntu(fluid))
        # print("eps: ", eps)

        # print(fluid.temperature)
        self.qwater = eps * fluid.capacity(fluid.temperature) * (t2phase - fluid.temperature)

        # print("Power: ", self.qwater)
        fluid.temperature += self.qwater / fluid.capacity(fluid.temperature)
        return fluid

    def solvebalance(self, inflow: Fluid, power):
        thermoparams = (self.uaf(inflow.vdot), inflow.capacity(inflow.temperature), inflow.temperature, power * 1000.0)
        # print(" Thermoparams: ", thermoparams)

        # print("Thermo Balance: ", self.balancehx(inflow.temperature, *thermoparams))

        x = fsolve(self.balancehx, inflow.temperature, args=thermoparams)  # inflow.temperature Startwert.
        # x = fsolve(KD27.balancehx, 34.0,
        #            args=(self.ua, inflow.capacity(), inflow.temperature, compressor.power_condenser()))

        print(x[0])
        self.twophase_temperature = x[0]

        return x

    def balancehx(self, t2phase, *thermoparams):  # ersetzt Condenser Balance
        # print("In balancehx: ", thermoparams)
        ua, capacity, temperature_water_in, q_condenser = thermoparams
        ntu = ua / capacity
        eps = 1 - math.exp(-ntu)
        # uaeff = eps * capacity / (1 - eps)
        # qwater = uaeff * (t2phase - temperature_water_in)
        qwater = eps * capacity * (t2phase - temperature_water_in)
        print("In balancehx: result: ", q_condenser - qwater, qwater, eps,  ntu, t2phase)
        return q_condenser - qwater

    def fithx2phase(self, inflow: Fluid, power, t2ph):
        thermoparams = (
            t2ph, inflow.capacity(inflow.temperature), inflow.temperature, power * 1000.0, self.number_of_plates,
            self.area_of_plate)
        # print("Thermo Balance: ", KD27.balancehx(34, *thermoparams))

        htc = fsolve(self.fitbalance, 1000, args=thermoparams)
        return htc

    def fitbalance(self, htc, *thermoparams):
        #  print("In fitbalance: ", thermoparams)
        t2phase, capacity, temperature_water_in, q_condenser, noplates, areaplate = thermoparams
        ua = htc * noplates * areaplate
        ntu = ua / capacity
        eps = 1 - math.exp(-ntu)
        # uaeff = eps * capacity / (1 - eps)
        # qwater = uaeff * (t2phase - temperature_water_in)
        qwater = eps * capacity * (t2phase - temperature_water_in)
        #  print("In fitbalance: result: ", q_condenser - qwater, qwater, eps, uaeff, ntu, htc, ua)
        return q_condenser - qwater

    def test_HX(self, inflow, power):   # power in kW
        return self.solvebalance(inflow, power)

    def loop_vdot_HX(self, inflow: Fluid, power):    # power in kW
        data = []
        v0 = inflow.vdot

        for iflow in [-0.5, -0.25, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.25, 0.5]:
            inflow.vdot = v0 * (1 + iflow)

            data.append([inflow.temperature, inflow.vdot, power, self.test_HX(inflow, power)[0],
                         inflow.heat(power*1000.0)])
        inflow.vdot = v0
        df = pd.DataFrame(data, columns=["Tsource", "vdot source", " Power", "T2phase", " Toutsource"])
        print(df)


class Evaporator(TwophaseHX):
    pass


class Evaporator_fitted(TwophaseHX):
    def __init__(self, no_plates, area_plate, htc, ua):
        super().__init__(no_plates, area_plate, htc)
        self.uacoeffs = ua



class Condenser(TwophaseHX):
    pass
    # def balance(self, t2phasex, compressor: Compressor, HX: HeatExchanger, fluid: Fluid):  # aktuell aus dem Objekt
    #     HX.calculate(fluid, t2phasex)
    #     return compressor.q_condenser() - HX.qwater


class Condenser_fitted(TwophaseHX):
    def __init__(self, no_plates, area_plate, htc, ua):
        super().__init__(no_plates, area_plate, htc)
        self.uacoeffs = ua


class Heater(HeatExchanger):
    def __init__(self, max_power, efficiency, costs, area, htc):
        self.efficiency = efficiency
        self.power_max = max_power
        self.costs_per_kwh = costs
        self.heat_transfer_coeff = htc
        super().__init__(1.0, area, htc)
        self.kW_per_kg_gas = 13.8


    def power(self, fraction):
        return self.power_max * fraction * 1000.0   # in Watt

    def costs(self, fraction):
        return self.power(fraction) * self.costs_per_kwh / self.efficiency

    def mass_exhaust(self, power):
        return 14.23 * 1.1 * power / self.kW_per_kg_gas / 3600.0  # Heizwert * 10% Luftueberschuss in kg/s

    def exhaust_gas_temperature(self, fraction):
        me = self.mass_exhaust(self.power(fraction) / 1000.0)
        hg = Burned_Air(0, me)
        return hg.heat(self.power(fraction))

    def power_water(self, hotgas, coldflow):
        return self.calculate(hotgas, coldflow, 1)

    def hot_flow(self, power, fraction):
        hg = Burned_Air(0,self.mass_exhaust(self.power_max * fraction))
        hg.temperature = hg.heat(self.power(fraction))
        return hg



class AirWaterunit(HeatExchanger):
    def __init__(self, fan: Fan, air: WetAir, area, htc):
        self.Fan = fan
        self.Air = air
        super().__init__(1.0, area, htc)

    def calculate_with_fan(self, air, water):
        qwater = self.calculate(air, water, 2)   # 2   crossflow
        return qwater + 0.9 * self.Fan.electric_power

    def set_air_temperature(self, temp):
        self.Air.temperature = temp

    def set_air_volume(self, vol, water: Fluid =""):
        if water == "":   # fixed speed set
            self.Air.vdot = vol
        else:  # optimise for balances heat capacity
            cpx = water.capacity(self.Air.temperature)
            rd = self.Air.density(self.Air.temperature)
            self.Air.vdot = cpx / self.Air.cp_t(self.Air.temperature) / rd
        return self.Air.vdot


class GasCooler(HeatExchanger):
    def __init__(self, pump, kA):
        self.Pump = pump
        self.ua = kA

    def power(self, heatwater_in: Fluid, t_hotgas):
        if heatwater_in.temperature > t_hotgas:
            heatwater_in.vdot = 0.0
            qheat = 0.0
        else:
            qheat = self.ua * (t_hotgas - heatwater_in.temperature) / 1000.0
        return qheat


class NoGasCooler(GasCooler):
    def __init__(self):
        pump = NoPump
        super(NoGasCooler, self).__init__(pump, 0.0)


