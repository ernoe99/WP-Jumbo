import math
import pandas as pd

from scipy.optimize import fsolve

from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.CoolProp import PropsSI, PhaseSI

from Fluids import Fluid, Heatwater, Brine, Air, WetAir, Burned_Air, HAir
from PumpsandFans import Fan

K0 = 273.15
p0 = 101325.0

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

class OnephaseHX_fitted(HeatExchanger):
    def __init__(self, no_plates, area_plate, htc, ua):
        super().__init__(no_plates, area_plate, htc)
        self.uacoeffs = ua

    def uaf(self, fluid: Fluid):

        if len(self.uacoeffs) == 3:
           uaft = (self.uacoeffs[0] + self.uacoeffs[1] * fluid.vdot + self.uacoeffs[2] * fluid.vdot**2) * \
                  self.area_of_plate * self.number_of_plates
        elif len(self.uacoeffs) == 8:
            uaft1 = (self.uacoeffs[1] + self.uacoeffs[2] * fluid.vdot + self.uacoeffs[3]* fluid.vdot**2)
            uaft2 = (self.uacoeffs[5] + self.uacoeffs[6] * fluid.vdot + self.uacoeffs[7] * fluid.vdot ** 2)

            uaft = (uaft1 + (uaft2 - uaft1) * (fluid.temperature - self.uacoeffs[0]) / (
                        self.uacoeffs[4] - self.uacoeffs[0])
                    ) * self.area_of_plate * self.number_of_plates
        else:
            uaft = self.ua

        return uaft

    def ntu(self, fluid: Fluid):
        return self.uaf(fluid) / fluid.capacity(fluid.temperature)

    def calculate(self, hotflow: Fluid, coldflow: Fluid, iflow):
        # ntu = max(self.ntu(hotflow), self.ntu(coldflow))
        print("NTU: ", self.ntu(hotflow), self.ntu(coldflow), self.ua, hotflow.capacity(hotflow.temperature),
              coldflow.capacity(coldflow.temperature))
        ntu = self.uaf(hotflow) / min(hotflow.capacity(hotflow.temperature), coldflow.capacity(coldflow.temperature))
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

    def uaf(self, fluid: Fluid):

        if len(self.uacoeffs) == 3:
            uaft = (self.uacoeffs[0] + self.uacoeffs[1] * fluid.vdot + self.uacoeffs[2] * fluid.vdot ** 2) * \
                   self.area_of_plate * self.number_of_plates
        elif len(self.uacoeffs) == 8:
            uaft1 = (self.uacoeffs[1] + self.uacoeffs[2] * fluid.vdot + self.uacoeffs[3] * fluid.vdot ** 2)
            uaft2 = (self.uacoeffs[5] + self.uacoeffs[6] * fluid.vdot + self.uacoeffs[7] * fluid.vdot ** 2)

            uaft = (uaft1 + (uaft2 - uaft1) * (fluid.temperature - self.uacoeffs[0]) / (
                    self.uacoeffs[4] - self.uacoeffs[0])
                    ) * self.area_of_plate * self.number_of_plates
        else:
            uaft = self.ua

        return uaft


    def ntu(self, fluid: Fluid):
        return self.uaf(fluid) / fluid.capacity(fluid.temperature)

    def calculate(self, fluid: Fluid, t2phase):
        eps = 1 - math.exp(-self.ntu(fluid))
        # print("eps: ", eps)

        # print(fluid.temperature)
        self.qwater = eps * fluid.capacity(fluid.temperature) * (t2phase - fluid.temperature)

        # print("Power: ", self.qwater)
        fluid.temperature += self.qwater / fluid.capacity(fluid.temperature)
        return fluid

    def solvebalance(self, inflow: Fluid, power):
        thermoparams = (self.uaf(inflow), inflow.capacity(inflow.temperature), inflow.temperature, power * 1000.0)
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


# Luftwärmepumpen Objekte

# The classes


class iwt:  # special for Refrigerant / refrigerant
    def __init__(self, kA):

        self.kA = kA
        self.power = -1.0e99

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

        self.power = -q_max * eps
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

        self.power = -q_max * eps  # negativ, da für Verdampfer in LuftKM_OEMWT.solvebalance verwendet

        return q_max * eps


class noiwt(iwt):
    def __init__(self):
        self.kA = 0.0
        self.power = 0.0

    def calculate(self, mk1, hk1, mk2, hk2, ref, pk1, pk2):
        return 0.0




noiwt = noiwt()
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


class AirEvaporator:  # special for air and evaporation
    def __init__(self, kA):

        self.kA = kA

    #    def calculate(self, vol, Tairin, Rin, Tevap, SH, href_in, refrigerant, href_out=0.0):  # Rin in %, vol in m3/h, T in Celsius

    def calculate(self, air: HAir, Tevap, refstate_I: list, refstate_O : list,  refrigerant):
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
            Tair_out = tair_out_dry - K0
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
        h1sat = PropsSI('H', 'T', Tevap + K0, 'Q', 1, refrigerant.name)

        href_in = refstate_I[3]

        mdot_ref = qfull / (hC_in - href_in)
        Vair_out = mdotair / rhoair_out

        return [mdot_ref * 3600, qfull, qact, qcond + qice, Vair_out * 3600.0, Tair_out, tair_out_dry - K0,
                hC_in]

    def calc(self, Tevap, air: HAir):
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

    def solvebalance(self, iwtact : iwt, air : HAir,  T2phase, power, refstate_1, refstate_3E, r32):
        thermoparams = (iwtact, air, power, refstate_1, refstate_3E, r32)
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

        iwt, air, power, refstate_1, refstate_3E, r32 = thermoparams
        # power = refstate_1[0] * refstate_1[3] - refstate_3[0] * refstate_3[3]
        # iwt = thermoparams[1]
        # evap = thermoparams[0]
        # air = themoparams[2]

        power_iwt = iwt.calc(refstate_3E, refstate_1, r32.name)
        power_evap = self.calc(t2phase, air)
        print("In balancehx Evapiwt; ", power_iwt, power_evap, t2phase, power - (power_evap + power_iwt))
        return power - (power_evap + power_iwt)

