import pyromat as pm
import pandas as pd

import numpy as np
import math
from scipy.optimize import fsolve


class Fluid:
    def __init__(self, temp, volflow):
        self.temperature = temp
        self.vdot = volflow
        self.cp = 4000.0   # only initializing

    def mdot(self, temp):
        return self.vdot * self.density(temp)

    def density(self, temp):
        return 1.0

    def enthalpy(self):
        return self.temperature * self.mdot(self.temperature) * self.cp

    def getvector(self):
        return [self.mdot(self.temperature), self.vdot, self.temperature, self.density(self.temperature), self.enthalpy(), self.cp]

    def capacity(self, temp):
        return self.mdot(temp) * self.cp

    def heat(self, qdot):
        return self.temperature + qdot / self.capacity(self.temperature)

class Heatwater(Fluid):
    def __init__(self, temp, volflow):
        super().__init__(temp, volflow)
        self.cp = 4200.0  # override cp for Heatwater

class Brine(Fluid):
    def __init__(self, temp, volflow):
        super().__init__(temp, volflow)
        self.cp = 3600.0  # override cp for Brine


class Air(Fluid):
    def __init__(self, temp, volflow):
        super().__init__(temp, volflow)
        self.cp = 1000.0

    def density(self, temp):
        return 101325.0 / (287 * (temp + 273.15))    # Achtung kg/m3

class WetAir(Air):
    def __init__(self, temp, volflow, rf):
        super().__init__(temp, volflow)
        self.cp = 1000.0
        self.rF = rf     # relative Feuchte

    def pds_t(self, temp):
        if temp < 0:
            pds_tx = 611.2 * (1 + temp / 148.6)**12.2
        else:
            pds_tx = math.exp(23.62-4065.0/(temp+236.27))
        return pds_tx

    def density_rF(self, temp, rF, p):
        # t in °C
        # rF relative Feuchte[0 - 1]
        # p in Pascal

        mDampf = 18.0  # g / mol
        mLuft = 28.96  # g / mol
        r = 8314.472  # J / kg / K

        phips = rF * self.pds_t(temp)
        return 1 / (r * (273.15 + temp)) * (mDampf * phips + mLuft * (p - phips))

    def cp_t(self, temp):
        return (1.0032+((1424/(100+max(-50,temp)))**1.4+1.94)**-2)*1000.0

    def feuchtkugel(self, temp, phi):
        return -5.809 + 0.058*phi*100 + 0.697*temp + 0.003*phi*100.0*temp

    def hv_t(self, temp):
        # % Verdampfungsenthalpie von Wasser
        # % t in Celsius
        # % hv_t in kJ / kg
        return 2500.9 / (1 + (0.2281 + (max(0, temp) / 4000) ** 1.15) * math.tan(temp / 238.1))

    def xs_t_p(self, temp, p):
        # % x bei 100 % Luftfeuchte
        # % t in Celsius
        # % p in Pa
        # % xs in kg / kg
        pds = self.pds_t(temp)
        if pds > p:
            xs_t_p = 20.0
        else:
            xs_t_p = 0.622 * pds / (p - pds)
        return xs_t_p

    def x_t_phi_p(self, temp, phi, p):
        pdx = self.pds_t(temp) * phi
        return 0.622 * pdx / (p - pdx)

    def taupunkt_t_phi(self, temp, phi):
        K1 = 6112.13
        K2 = 17.5043
        K3 = 241.2

        return K3 * ((K2 * temp) / (K3 + temp) + math.log(phi)) / (K2 * K3 / (K3 + temp) - math.log(phi));

    def t_h_x_p(self, h, x, p):   # Umkehrfunktion zu h_t_x_p
        return (h - x * self.hv_t(0.0)) / (1.006 + 1.86*x)

    def phi_t_x_p(self, temp, xact, p):
        pds = self.pds_t(temp)
        pd = p * x / (0.622 + xact)
        return min(pd / pds, 1.0)

    def h_t_x_p(self, temp, xact, p):
        xs = self.xs_t_p(temp, p)
        cpx = 1006.0
        if xact > xs:
            h = temp * cpx + xs * (hv_t(0) + 1.86 * temp)
            xliquid = xact - xs
        else:
            h = temp * cpx + x * (hv_t(0) + 1.86 * temp)
            xliquid = 0
        return [h, xliquid]

    def condens(self, t_target, x, p):
        xs = self.xs_t_p(t_target, p)
        return max((x - xs), 0.0)

class Sink:  # This is the building
    def __init__(self, heatload):
        self.q_heat = heatload

    def calculate(self, sink_in : Fluid, sink_out : Fluid):
        sink_out.temperature = (sink_in.temperature * sink_in.cp - self.q_heat) / sink_in.cp
        return (sink_in.temperature * sink_in.cp - self.q_heat) / sink_in.cp

    def rl_temperature(self, sink_in: Fluid):
        # print("rl_temperature", sink_in.temperature, sink_in.capacity(), self.q_heat)
        return (sink_in.temperature * sink_in.capacity(sink_in.temperature) - self.q_heat) / sink_in.capacity(sink_in.temperature)

class Pump:
    def __init__(self, max_flow):
        self.max_flow = max_flow

    def flow(self, fraction):
        return self.max_flow * fraction


class EXV:
    def __init__(self, superheat):
        self.superheat = superheat

    def calculate(self, sensor: Fluid):
        pass  # TODO


class Compressor:
    def __init__(self, hsuc, hdisch, hliquid, mflow):
        self.enthalpy_liquid_subcooled = hliquid
        self.enthalpy_suction = hsuc
        self.enthalpy_discharge = hdisch
        self.massflow_discharge = mflow

    def q_condenser(self):
        return (self.enthalpy_discharge - self.enthalpy_liquid_subcooled) * self.massflow_discharge

    def calculate(self, input):
        power, ssh, sdt, sst, llt, economizer, tapp = input

class TurboCore(Compressor):
    def __init__(self, fname : str):
        # optimal data w/o Economizer
        fcsv = fname + '-noEcon_0Power.csv'
        pddata = pd.read_csv(fcsv, header=9, sep=';', nrows=527)
        linedat = pddata.to_numpy(copy=True)
        self.zeropower_noecon_data = linedat.reshape(17, 31, 39)

        fcsv = fname + '-Econ_0Power.csv'
        pddata = pd.read_csv(fcsv, header=9, sep=';', nrows=527)
        linedat = pddata.to_numpy(copy=True)
        self.zeropower_econ_data = linedat.reshape(17, 31, 39)

        fcsv = fname + '-noEcon_FromPower.csv'
        pddata = pd.read_csv(fcsv, header=9, sep=';', nrows=2635)
        linedat = pddata.to_numpy(copy=True)
        self.frompower_noecon_data = linedat.reshape(5, 17, 31, 39)

        fcsv = fname + '-Econ_FromPower.csv'
        pddata = pd.read_csv(fcsv, header=9, sep=';', nrows=2635)
        linedat = pddata.to_numpy(copy=True)
        self.frompower_econ_data = linedat.reshape(5, 17, 31, 39)

        self.actual_values = np.empty(38)
        self.economizer = 0

    def calculate_0power_econ(self, t_suction, t_condensation):

        id = int((t_suction + 18) / 3)
        jd = int((t_condensation + 5) / 2.5)
        vec0 = np.empty(38)

        if -17.99 <= t_suction <= 30 and -5 <= t_condensation <= 69.99:
            rs = self.zeropower_econ_data

            if rs[id, jd, 14] != 0 and rs[id + 1, jd, 14] != 0 and rs[id, jd + 1, 14] != 0 and rs[
                id + 1, jd + 1, 14] != 0:

                for i in range(1, 39, 1):  # linear interpolation
                    w11 = rs[id, jd, i] + (rs[id + 1, jd, i] - rs[id, jd, i]) / 3.0 * (t_suction - rs[id, jd, 1])
                    w12 = rs[id, jd + 1, i] + (rs[id + 1, jd + 1, i] - rs[id, jd + 1, i]) / 3.0 * (
                            t_suction - rs[id, jd + 1, 1])
                    w23 = w11 + (w12 - w11) / 2.5 * (t_condensation - rs[id, jd, 2])
                    vec0[i - 1] = w23
                if vec0[10] != 0:
                    self.actual_values = vec0
                    valid = 1
            else:
                #   massflow zero means no valid result from compressor in this point
                valid = 0
        else:
            valid = 0  # temperatures out of range
        return valid

    def calculate_0power_noecon(self, t_suction, t_condensation):

        id = int((t_suction + 18) / 3)
        jd = int((t_condensation + 5) / 2.5)
        vec0 = np.empty(38)

        if -17.99 <= t_suction <= 30 and -5 <= t_condensation <= 69.99:

            rs = self.zeropower_noecon_data

            if rs[id, jd, 14] != 0 and rs[id + 1, jd, 14] != 0 and rs[id, jd + 1, 14] != 0 and rs[
                id + 1, jd + 1, 14] != 0:
                for i in range(1, 39, 1):  # linear interpolation
                    w11 = rs[id, jd, i] + (rs[id + 1, jd, i] - rs[id, jd, i]) / 3.0 * (t_suction - rs[id, jd, 1])
                    w12 = rs[id, jd + 1, i] + (rs[id + 1, jd + 1, i] - rs[id, jd + 1, i]) / 3.0 * (
                            t_suction - rs[id, jd + 1, 1])
                    w23 = w11 + (w12 - w11) / 2.5 * (t_condensation - rs[id, jd, 2])
                    vec0[i - 1] = w23
                self.actual_values = vec0
                valid = 1
            else:
                #   massflow zero means no valid result from compressor in this point
                valid = 0
        else:
            valid = 0  # temperatures out of range
        return valid

    def calculate_0power(self, t_suction, t_condensation):
        valid = self.calculate_0Power_Econ(t_suction, t_condensation)
        if valid == 0:
            valid = self.calculate_0Power_noEcon(t_suction, t_condensation)
        return valid  # Results in self.actual_values_0Power

    def get_values(self, t_suction, t_condensation, rs):
        id = int((t_suction + 18) / 3)
        jd = int((t_condensation + 5) / 2.5)
        vec0 = np.empty(38)

        if -17.99 <= t_suction <= 30 and -5 <= t_condensation <= 69.99:

            if rs[id, jd, 14] != 0 and rs[id + 1, jd, 14] != 0 and rs[id, jd + 1, 14] != 0 and rs[id + 1, jd + 1, 14] != 0:
                for i in range(1, 39, 1):  # linear interpolation
                    w11 = rs[id, jd, i] + (rs[id + 1, jd, i] - rs[id, jd, i]) / 3.0 * (t_suction - rs[id, jd, 1])
                    w12 = rs[id, jd + 1, i] + (rs[id + 1, jd + 1, i] - rs[id, jd + 1, i]) / 3.0 * (
                            t_suction - rs[id, jd + 1, 1])
                    w23 = w11 + (w12 - w11) / 2.5 * (t_condensation - rs[id, jd, 2])
                    vec0[i - 1] = w23
                self.actual_values = vec0
                valid = 1
            else:
                #   massflow zero means no valid result from compressor in this point
                valid = 0
        else:
            valid = 0  # temperatures out of range
        return valid

    def calculate_zeropower(self, t_suction, t_condensation):

        valid = self.get_values(t_suction, t_condensation, self.zeropower_econ_data)

        if valid == 0:
            valid = self.get_values(t_suction, t_condensation, self.zeropower_noecon_data)
            self.economizer = 0
        else:
            self.economizer = 1
        return valid


    def calculate_frompower(self, frompower, t_suction, t_condensation):

        id1 = int(frompower / 25.0)
        valid1 = self.get_values(t_suction, t_condensation, self.frompower_econ_data[id1])
        vec1 = self.actual_values
        valid2 = self.get_values(t_suction, t_condensation, self.frompower_econ_data[min(id1 + 1, 4)])
        vec2 = self.actual_values
        self.economizer = 1

        if valid1 == 0 or valid2 == 0:  # probiere ohne Economizer
            valid1 = self.get_values(t_suction, t_condensation, self.frompower_noecon_data[id1])
            vec1 = self.actual_values
            valid2 = self.get_values(t_suction, t_condensation, self.frompower_noecon_data[min(id1 + 1, 4)])
            vec2 = self.actual_values
            self.economizer = 0

        if valid1 == 1 and valid2 == 1:
            ref = [0.0, 25.0, 50.0, 75.0, 100.0]
            for i in range(0, 38, 1):  # linear interpolation achtung 1 niedriger da, CID schon geschnitten
                vec1[i] = vec1[i] + (vec2[i] - vec1[i]) / 25.0 * (frompower - ref[id1])
            self.actual_values = vec1
            valid = 1
        else:
            valid = 0
        return valid

    def power_evaporator(self):
        return self.actual_values[17]

    def power_electrical(self):
        return self.actual_values[16]

    def power_condenser(self):
        return (self.actual_values[10] - self.actual_values[9]) * (self.actual_values[13] + self.actual_values[32])

    def cop_cooling(self):
        return self.actual_values[15]

    def cop(self):
        return self.actual_values[15] + 1.0

    def massflow_evaporator(self):
        return self.actual_values[13]

    def massflow_condenser(self):
        return self.actual_values[13] + self.actual_values[32]

    def pressure_ratio(self):
        return self.actual_values[25]

    def temperature_discharge(self):
        return self.actual_values[24]

    def temperature_suction(self):
        return self.actual_values[0] + self.actual_values[20]

    def pressure_suction(self):
        return self.actual_values[22] / 100.0  # Umrechnung in bar

    def pressure_discharge(self):
        return self.actual_values[23] / 100.0  # Umrechnung in bar

    def power_economizer(self):
        return self.actual_values[31]

    def pressure_economizer_interstage(self):
        return self.actual_values[30] / 100.0  # Umrechnung in bar

    def temperature_economizer_interstage(self):
        return self.actual_values[29]

    def massflow_economizer(self):
        return self.actual_values[32]


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

    def calculate(self, hotflow : Fluid, coldflow : Fluid, iflow):
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

        if iflow == 1:   # Gegenstrom
            if cr > 0.999:  # Sonderfall
                eps = ntu / (1 + ntu)
            else:
                eps = (1 - math.exp(-ntu * (1 - cr))) / (1 - cr * math.exp(-ntu * (1 - cr)))
        if iflow == 2:    # Kreuzstrom
            #  eps = 1 - math.exp(qmax**0.22/cr * (math.exp(-cr * qmax**0.78) - 1))
            eps = 1 - math.exp((math.exp(-cr * ntu**0.78) - 1) / (cr * ntu**-0.22))

            print((cr * ntu**-0.22))
            print(math.exp(-cr * ntu**0.78) - 1)
            print(math.exp(-cr * ntu**0.78))
            print(-cr)
            print(ntu**0.78)
            print(math.exp((math.exp(-cr * ntu**0.78) - 1) / (cr * ntu**-0.22)))

        print("EPS: ", eps, q_max * eps)
        return q_max * eps



    def solvebalance(self, inflow: Fluid, power):
        thermoparams = (self.ua, inflow.capacity(inflow.temperature), inflow.temperature, power*1000.0)
        # print(" Thermoparams: ", thermoparams)

        print("Thermo Balance: ", self.balancehx(inflow.temperature, *thermoparams))

        x = fsolve(self.balancehx, inflow.temperature, args=thermoparams)  # inflow.temperature Startwert.
        # x = fsolve(KD27.balancehx, 34.0,
        #            args=(self.ua, inflow.capacity(), inflow.temperature, compressor.power_condenser()))

        print(x[0])
        self.twophase_temperature = x[0]

        return x

    def balancehx(self, t2phase, *thermoparams):  # ersetzt Condenser Balance
        print("In balancehx: ", thermoparams)
        ua, capacity, temperature_water_in, q_condenser = thermoparams
        ntu = ua / capacity
        eps = 1 - math.exp(-ntu)
        uaeff = eps * capacity / (1 - eps)
        qwater = uaeff * (t2phase - temperature_water_in)
        print("In balancehx: result: ", q_condenser - qwater, qwater, eps, uaeff, ntu, t2phase)
        return q_condenser - qwater

    def fithx2phase(self, inflow: Fluid, power, t2ph):
        thermoparams = (t2ph, inflow.capacity(inflow.temperature), inflow.temperature, power * 1000.0, self.number_of_plates,
                        self.area_of_plate)
        # print("Thermo Balance: ", KD27.balancehx(34, *thermoparams))

        htc = fsolve(self.fitbalance, 1000, args=thermoparams)
        return htc

    def fitbalance(self, htc, *thermoparams):
        print("In fitbalance: ", thermoparams)
        t2phase, capacity, temperature_water_in, q_condenser, noplates, areaplate = thermoparams
        ua = htc * noplates * areaplate
        ntu = ua / capacity
        eps = 1 - math.exp(-ntu)
        uaeff = eps * capacity / (1 - eps)
        qwater = uaeff * (t2phase - temperature_water_in)
        print("In fitbalance: result: ", q_condenser - qwater, qwater, eps, uaeff, ntu, htc, ua)
        return q_condenser - qwater


class TwophaseHX(HeatExchanger):
    def calculate(self, fluid: Fluid, t2phase):
        eps = 1 - math.exp(-self.ntu(fluid))
        print("eps: ", eps)
        uaeff = eps * fluid.capacity(fluid.temperature) / (1 - eps)
        print("uaeff: ", uaeff)
        print(t2phase)

        print(fluid.temperature)
        self.qwater = uaeff * (t2phase - fluid.temperature)
        print("Power: ", self.qwater)
        fluid.temperature += self.qwater / fluid.capacity(fluid.temperature)
        return fluid


class Evaporator(TwophaseHX):
    pass


class Condenser(TwophaseHX):
    pass
    # def balance(self, t2phasex, compressor: Compressor, HX: HeatExchanger, fluid: Fluid):  # aktuell aus dem Objekt
    #     HX.calculate(fluid, t2phasex)
    #     return compressor.q_condenser() - HX.qwater


class Heater:
    def __init__(self, max_power, efficiency, costs, htc):
        self.efficiency = efficiency
        self.power_max = max_power
        self.costs_per_kwh = costs
        self.heat_transfer_coeff = htc

    def power(self, fraction):
        return self.power_max * fraction

    def costs(self, fraction):
        return self.power(fraction) * self.costs_per_kwh


class HeatPump:
    def __init__(self, compressor: Compressor, condenser: Condenser, evaporator: Evaporator, source_in: Fluid,
                 source_out: Fluid, sink_in: Fluid, sink_out: Fluid):
        self.components = (sink_in, compressor, condenser, sink_out, source_in, evaporator, source_out)

class SimpleHeatPump:
    def __init__(self, compressor: TurboCore, condenser: Condenser, evaporator: Evaporator, source_in: Fluid,
                 source_out: Fluid, sink_in: Fluid, sink_out: Fluid):
        self.target_temperature = sink_out.temperature
        self.TurboCore = compressor
        self.condenser = condenser
        self.evaporator = evaporator
        self.source_in  = source_in
        self.source_out = source_out
        self.sink_in = sink_in
        self.sink_out = sink_out
        print("SimpleHeatPump initialized", compressor)

    def calculate(self, t_target, frompower, t_suction, t_condensation):
        valid = self.TurboCore.calculate_frompower(frompower, t_suction, t_condensation)
        if valid != 1:
            print("Failure in Compressor valid = ", valid)
            quit("Failure")
        t_condensation_new = self.condenser.solvebalance(self.sink_in, self.TurboCore.power_condenser())
        print("TurboCore Condenser Power: ", self.TurboCore.power_condenser(), t_condensation_new)
        self.sink_out.temperature = self.sink_in.heat(self.TurboCore.power_condenser()*1000.0)
        print("Sink out Temperature: ", self.sink_out.temperature)
        t_suction_new = self.evaporator.solvebalance(self.source_in, -self.TurboCore.power_evaporator())
        print("TurboCore Evaporator Power: ", self.TurboCore.power_evaporator(), t_suction_new)
        self.source_out.temperature = self.source_in.heat(-TurboCore.power_evaporator(self.TurboCore)*1000.0)
        print("Temperature Source_out: ", self.source_out.temperature)
        return [frompower, t_suction_new, t_condensation_new, t_target - self.sink_out.temperature]

class HeatPumpEconomizer(HeatPump):
    def __init__(self, compressor: Compressor, condenser: Condenser, evaporator: Evaporator, source_in: Fluid,
                 source_out: Fluid, sink_in: Fluid, sink_out: Fluid, economizer: Evaporator):
        self.components = (sink_in, economizer, compressor, condenser, sink_out, source_in, evaporator, source_out)

class TwoStageHeatPump(HeatPump):
    def __init__(self, hp1:HeatPump, hp2: HeatPump):
        self.components = (hp1.components + hp2.components)

# Umwelt

aussentemperatur = 2.0

# Beginn Aufbau System

vorlauftemperatur = 35.0
vdot_init = 18.0
qload_sink = 250000.0

heizkreispumpe_maxflow = 25.0
heizkreispumpe_actor = 18/25.0     # control

#Heatpump
condenser_no_plates = 120.0
condenser_area_plate = 0.15
condenser_heat_transfer = 640  # TODO Liste mit Komponenten
evaporator_no_plates = 240.0
evaporator_area_plate = 0.25
evaporator__heat_transfer = 440  # TODO Liste mit Komponenten

#compressor

sat_discharge_temperature = vorlauftemperatur + 2
sat_suction_temperature = aussentemperatur - 8
economizer_approach_temperature = 0.0
superheat = 4.0


#intialisiere Heizkreis

sf = Fluid(35, 18)
print(sf.temperature)

Heizkreispumpe = Pump(heizkreispumpe_maxflow)
print(Heizkreispumpe.flow(heizkreispumpe_actor))

Vorlauf_Sink = Brine(vorlauftemperatur, Heizkreispumpe.flow(heizkreispumpe_actor))  # Interface Vorlauf

Gebaude = Sink(qload_sink)
print(Gebaude.rl_temperature(Vorlauf_Sink))

Ruecklauf_Sink = Brine(Gebaude.rl_temperature(Vorlauf_Sink), Heizkreispumpe.flow(heizkreispumpe_actor))

#initialisiere Heatpump

KD27 = Condenser(condenser_no_plates, condenser_area_plate, condenser_heat_transfer)

TGH285 = TurboCore('TGH285')

af = Fluid(0, 25.0)

EV27 = Evaporator(evaporator_no_plates, evaporator_area_plate, evaporator__heat_transfer)

# Initialisiere Ausseneinheit



# Berechnung
frompower = 50.0
valid = TGH285.calculate_frompower(frompower, sat_suction_temperature, sat_discharge_temperature)
print(valid)

print(TGH285.actual_values)
print(Ruecklauf_Sink)

ka = KD27.ua

for i in range(1, 10):

    print("IN LOOP", i, ka * i)
    KD27.ua = ka * i
    x = KD27.solvebalance(Ruecklauf_Sink, TGH285.power_condenser())
    print("X: ", x)
    print("KD27 two phase temperature: ", KD27.twophase_temperature)
    new_vorlauf_temperature = Ruecklauf_Sink.heat(TGH285.power_condenser() * 1000.0)
    print("Vorlauf neu: ", new_vorlauf_temperature)


y = EV27.solvebalance(af, -TGH285.power_evaporator())

print(KD27.ua)
print("NTU_Heatwater: ", KD27.ntu(Ruecklauf_Sink))

t2phase = 34.0
print(Ruecklauf_Sink.temperature)

watercondout = KD27.calculate(Ruecklauf_Sink, t2phase)

print(watercondout.temperature)
print(watercondout.getvector())

qcondenser_comp = (238367.0 - 46636.0) * (1.023 + 0.162)
c1 = Compressor(hdisch=238367.0, hsuc= 100000.0, hliquid=46636.0, mflow=1.023 + 0.162)

print(c1.massflow_discharge)

print("Qcondenser_comp: ", qcondenser_comp)
print("Qwater: ", KD27.qwater)


thermoparams = (KD27.ua, Ruecklauf_Sink.capacity(Ruecklauf_Sink.temperature), Ruecklauf_Sink.temperature, qcondenser_comp)
print(thermoparams)

# print("Thermo Balance: ", KD27.balancehx(34, *thermoparams))

# x = fsolve(KD27.balancehx, 34.0, args=thermoparams)
# x = fsolve(KD27.balancehx, 34.0, args=(KD27.ua, Ruecklauf_Sink.capacity(), Ruecklauf_Sink.temperature, qcondenser_comp))

x = KD27.solvebalance(Ruecklauf_Sink, TGH285.power_condenser())

# print(x[0])

print("KD27 two phase temperature: ", KD27.twophase_temperature)

tgh285 = TurboCore('TGH285')

print(tgh285.calculate_zeropower(-10, 25.1))
print(tgh285.actual_values)

print(tgh285.calculate_frompower(25.3, -5, 27.8))
print(tgh285.actual_values)

print(tgh285.power_evaporator())
print(tgh285.power_condenser())

print(tgh285.cop_cooling())
print(tgh285.power_electrical() + tgh285.power_evaporator())

res = np.empty([6, 10])

for i in range(0, 10, 1):
    valid = tgh285.calculate_frompower(i * 10.0, -5, 27.8)
    print(tgh285.economizer)
    res[0, i] = tgh285.power_evaporator()
    res[1, i] = tgh285.power_condenser()
    res[2, i] = tgh285.power_electrical()
    res[3, i] = tgh285.power_economizer()
    res[4, i] = tgh285.temperature_discharge()
    res[5, i] = tgh285.cop()


print(res.transpose())
frompower = 100.0
valid = TGH285.calculate_frompower(frompower, -6, 35)
print(valid)
print("Condenser Power 300.0: ", TGH285.power_condenser())


print(TGH285.actual_values)

AL_KB210_134AH_F = Condenser(134, 0.2485, 937.0)  # fitted

xf = Brine(28, 51652.0/3600)

ufit = AL_KB210_134AH_F.fithx2phase(xf, 300.0, 35.0)

print("UFIT: ", ufit)

print("NEW Temperature: ", xf.heat(300000.0))

hf = Brine(0.0, 61134/3600.0)
AL_ACH1000DQ_306AH = Evaporator(306, 0.4183, 233.0)  # fitted
hfit = AL_ACH1000DQ_306AH.fithx2phase(hf, -250.0, -6.5)

print("HFIT: ", hfit)
print("NEW Temperature: ", hf.heat(-250000.0))


air = Air(400, 2.8)
print("Air.massflow", air.mdot, air.density(air.temperature))
win = Fluid(40, 18)

hx = HeatExchanger(3, 1.5, 1000.0)
hx.calculate(air, win, 2)

hx.calculate(air, win, 1)

vol11 = 275000.0/3600.0/(12-7)
print("vol11: ", vol11)

f11 = Brine(12, vol11)
AC500EQ = Evaporator(306, 0.4183, 250)
hfit = AC500EQ.fithx2phase(f11, -275.0, 5.68)
print(" AC500EQ.ua: ", AC500EQ.ua)

print("HFIT: ", hfit)
print("NEW Temperature: ", f11.heat(-275000.0))

f12 = Fluid(26.0, 17.5)
CCF2880315_6P = Condenser(134, 0.2485, 937.0)
print("CCF2880315_6P.ua: ", CCF2880315_6P.ua)


print(f12.getvector())
f12.vdot = 100.0
print(f12.getvector())

#Initialize Heatpump

hrl = Heatwater(30, 15.0)
request_L = Heatwater(35.0, hrl.vdot)   # REquest from system 25 l/s mit 35 Grad
print("Requested Power Condenser: ", hrl.vdot, hrl.cp, hrl.vdot * hrl.cp * (request_L.temperature - hrl.temperature))
source_in_B0 = Brine(0, 25.0)       # Standard source for Thermalia 0 Grad
source_out = Brine(-4, source_in_B0.vdot)

t_suc = -6.0
t_cond = 36.0


Thermalia_275 = SimpleHeatPump(compressor=TGH285, condenser=AL_KB210_134AH_F, evaporator=AC500EQ, sink_in= hrl,
                               sink_out= request_L, source_in= source_in_B0, source_out=source_out)

val_iter = Thermalia_275.calculate(t_target=35.0, frompower=50.0, t_suction=t_suc, t_condensation= t_cond)
print("Iteration Values: ", val_iter)

# print(KD27.ua, KD27.number_of_plates, KD27.heat_transfer_coeff, KD27.area_of_plate)
# print(AC500EQ.ua, AC500EQ.number_of_plates, AC500EQ.heat_transfer_coeff, AC500EQ.area_of_plate)
# print(AL_KB210_134AH_F.ua, AL_KB210_134AH_F.number_of_plates, AL_KB210_134AH_F.heat_transfer_coeff,
#       AL_KB210_134AH_F.area_of_plate)
