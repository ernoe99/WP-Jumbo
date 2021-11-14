from TurboCor import Compressor, TurboCor
from HeatExchanger import HeatExchanger, Heater, Condenser, Evaporator_fitted, Evaporator, AirWaterunit, GasCooler, \
    Condenser_fitted
from Fluids import Fluid
import pandas as pd
import openpyxl


class HeatPump:
    def __init__(self, compressor: Compressor, condenser: Condenser, evaporator: Evaporator, source_in: Fluid,
                 source_out: Fluid, sink_in: Fluid, sink_out: Fluid):
        self.components = (sink_in, compressor, condenser, sink_out, source_in, evaporator, source_out)


class SimpleHeatPump:
    def __init__(self, compressor: TurboCor, condenser: Condenser_fitted, evaporator: Evaporator_fitted,
                 source_in: Fluid, source_out: Fluid, sink_in: Fluid, sink_out: Fluid):
        self.target_temperature = sink_out.temperature
        self.TurboCore = compressor
        self.condenser = condenser
        self.evaporator = evaporator
        self.source_in = source_in
        self.source_out = source_out
        self.sink_in = sink_in
        self.sink_out = sink_out
        print("SimpleHeatPump initialized", compressor)

    def p_condenser(self):
        return self.TurboCore.power_condenser()

    def p_evaporator(self):
        return self.TurboCore.power_evaporator()

    def calculate(self, frompower, t_suction, t_condensation):

        valid = self.TurboCore.calculate_frompower(frompower, t_suction, t_condensation)
        if valid != 1:
            print("Failure in Compressor valid = ", valid)
            quit("Failure")
        t_condensation_new = self.condenser.solvebalance(self.sink_in, self.p_condenser())[0]
        print("TurboCore Condenser Power: ", self.p_condenser(), t_condensation_new)
        self.sink_out.temperature = self.sink_in.heat(self.p_condenser() * 1000.0)
        print("Sink out Temperature: ", self.sink_out.temperature)
        t_suction_new = self.evaporator.solvebalance(self.source_in, -self.p_evaporator())[0]
        print("TurboCore Evaporator Power: ", self.p_evaporator(), t_suction_new)
        self.source_out.temperature = self.source_in.heat(-self.p_evaporator() * 1000.0)
        print("Temperature Source_out: ", self.source_out.temperature)
        return [frompower, t_suction_new, t_condensation_new, self.sink_out.temperature]

    def get_electrical_power(self):
        return self.TurboCore.power_electrical()

    def calculate_from_power(self, percent_power):
        # Simply give results or return 0
        iflop = 0

        t_suction = self.source_in.temperature - 3.0
        t_condensation = self.sink_in.temperature + 1.0
        t_suction_new = 999.0
        t_condensation_new = 999.0
        valid = self.TurboCore.calculate_frompower(percent_power, t_suction, t_condensation)
        if valid != 1:
            print("Failure in Compressor valid = ", valid)
            return valid

        while abs(t_suction - t_suction_new) and abs(t_condensation - t_condensation_new) > 0.1 and iflop < 99:
            print(" =====", abs(t_suction - t_suction_new), abs(t_condensation - t_condensation_new), iflop)
            if iflop != 0:
                t_suction = t_suction * 0.2 + t_suction_new * 0.8
                t_condensation = t_condensation * 0.7 + t_condensation_new * 0.3
            iflop += 1
            if iflop == 99:
                print("******", iflop, t_suction, t_condensation)

            [percent_power, t_suction_new, t_condensation_new, tsinkout] = \
                self.calculate(percent_power, t_suction, t_condensation)

            valid = self.TurboCore.calculate_frompower(percent_power, t_suction_new, t_condensation_new)
            if valid != 1:
                print("Failure in Compressor valid = ", valid)
                return valid

        return [percent_power, t_suction, t_condensation,
                self.sink_out.temperature, self.p_condenser()]  # nicht _new,da sonst Kompressordaten alt

    def calculate_from_flow(self, percent_power, t_target):
        # Simply give results or return 0

        isolve = 0
        iter = 0
        t_suction = self.source_in.temperature - 1.0
        t_condensation = self.sink_in.temperature + 2.0
        t_suction_new = 999.0
        t_condensation_new = 999.0
        valid = self.TurboCore.calculate_frompower(percent_power, t_suction, t_condensation)
        if valid != 1:
            print("Failure in Compressor valid = ", valid)
            return valid

        while isolve == 0:
            valid = self.TurboCore.calculate_frompower(percent_power, t_suction, t_condensation)
            if valid != 1:
                print("Failure in Compressor valid = ", valid)
                return valid

            [percent_power, t_suction_new, t_condensation_new, tsinkout] = self.calculate(percent_power, t_suction,
                                                                                          t_condensation)

            # t_condensation_new = self.condenser.solvebalance(self.sink_in, self.TurboCore.power_condenser())[
            #     0]  # solvebalance return x[0]
            # print("TurboCore Condenser Power: ", self.TurboCore.power_condenser(), t_condensation_new)
            # self.sink_out.temperature = self.sink_in.heat(self.TurboCore.power_condenser() * 1000.0)
            # print("Sink out Temperature: ", self.sink_out.temperature)
            # t_suction_new = self.evaporator.solvebalance(self.source_in, -self.TurboCore.power_evaporator())[0]
            # print("TurboCore Evaporator Power: ", self.TurboCore.power_evaporator(), t_suction_new)
            # self.source_out.temperature = self.source_in.heat(-TurboCor.power_evaporator(self.TurboCore) * 1000.0)
            # print("Temperature Source_out: ", self.source_out.temperature)

            if max(abs(t_suction - t_suction_new), abs(t_condensation - t_condensation_new), abs(t_target - self.sink_out.temperature)) > 0.05:
                isolve = 0
                print(" =====", abs(t_suction - t_suction_new), abs(t_condensation - t_condensation_new), abs(t_target - self.sink_out.temperature), iter)
                t_suction = t_suction_new
                t_condensation = t_condensation_new
                self.sink_in.vdot = self.sink_in.vdot * (self.sink_out.temperature - self.sink_in.temperature) \
                                    / (t_target - self.sink_in.temperature)
                self.sink_out.vdot = self.sink_in.vdot
                iter += 1
                print("NEW Volume: ", self.sink_in.vdot, t_target, self.sink_out.temperature)
                if iter > 99:
                    isolve = 1
                    with open('Errorlog.txt', 'a') as wfile:
                        wfile.write("*********************" + str(iter) + "  ABBRUCH 100 Iterationen" + "\n")
                        wfile.write(f" =====  {(abs(t_suction - t_suction_new))} ====  "
                                    f"{(abs(t_condensation - t_condensation_new))} "
                                    f"{(abs(t_target - self.sink_out.temperature))}  \n")
                        wfile.write(f"{self.sink_in.vdot}  === {t_target}   ** {self.source_in.vdot} \n")
            else:
                isolve = 1

        return [percent_power, t_suction, t_condensation,
                self.sink_out.temperature, self.TurboCore.power_condenser()]  # nicht _new,da sonst Kompressordaten alt

    def write_log(self, x):   # x is result from from_flow or from_power
        filename = "SimpleHeatPump" + ".log"

        if x == 0:
            with open(filename, 'a') as log:
                log.write('Compressor Failure: t_suction t_condensation from_power' + str(x) + "\n \n")
        else:
            with open(filename, 'a') as log:
                comp = self.TurboCore
                log.write("Power: " + str(x) + "\n")
                log.write("Power Data from Compressor:{0} {1} {2} {3} {4}\n".format(str(comp.power_electrical()),
                                                                                    str(comp.power_evaporator()),
                                                                                    str(comp.power_condenser()),
                                                                                    str(comp.power_economizer()),
                                                                                    str(comp.power_condenser() /
                                                                                        comp.power_electrical())))
                x = [self.sink_out.getvector(), self.sink_in.temperature,
                     self.sink_out.enthalpy() - self.sink_in.enthalpy()]
                log.write("Sink: " + str(x) + "\n")

                x = [self.source_out.getvector(), self.source_in.temperature,
                     self.source_out.enthalpy() - self.source_in.enthalpy()]
                log.write("Source: " + str(x) + "\n" + "\n")

    def operating_field_single(self, mode, target, filename):
        #  mode = 1 From_Flow    target = Vorlauftemperatur
        #  mode = 2 From Power   target = power in percent

        v0sink = self.sink_in.vdot
        v0source = self.source_in.vdot

        if mode == 1:
            request_VL = target
        if mode == 2:
            requested_Power = target
        
        data = []
        nosolution = []

        for ivol_sink in [-0.5, -0.25, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5]:
            for ivol_source in [-0.5, -0.25, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5]:
                for ipower in [0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:

                    self.sink_in.vdot = v0sink * (ivol_sink + 1)
                    self.sink_out.vdot = self.sink_in.vdot
                    self.source_in.vdot = v0source * (ivol_source + 1)
                    self.source_out.vdot = self.source_in.vdot

                    if mode == 1:
                        x = self.calculate_from_flow(ipower, request_VL)
                    else:
                        x = self.calculate_from_power(ipower, requested_Power)

                    print(v0sink * (ivol_sink + 1), v0source * (ivol_source + 1), ipower)

                    if x != 0:
                        data.append([self.source_in.vdot, self.source_in.temperature, self.source_out.temperature,
                                     self.sink_in.vdot, self.sink_in.temperature, self.sink_out.temperature, ipower,
                                     self.TurboCore.power_condenser(),
                                     self.TurboCore.power_condenser() / self.TurboCore.power_electrical(),
                                     self.TurboCore.power_evaporator(), self.TurboCore.power_electrical(),
                                     self.TurboCore.temperature_suction(), self.TurboCore.temperature_discharge_saturated(),
                                     self.TurboCore.actual_values[0]
                                     ])
                        print("Power: ", self.TurboCore.power_condenser(),
                              self.TurboCore.power_condenser() / self.TurboCore.power_electrical())
                        self.write_log(x)
                    else:
                        nosolution.append([self.source_in.vdot, self.source_in.temperature, self.source_out.temperature,
                                           self.sink_in.vdot, self.sink_in.temperature, self.sink_out.temperature, ipower])

        dfout = pd.DataFrame(data,
                             columns=["source Vdot", "Source Tin", "Source_Tout", "Sink Vdot", "Sink Tin", "Sink Tout",
                                      "Power Set", "Power Cond", "COPH", "Power Evap", "Power Electric", "Tsuction",
                                      " T Disch sat", "Actual Values"])
        dfout = dfout.sort_values(by='Power Cond', ascending=False)
        dfeff = dfout.sort_values(by='COPH', ascending=False)

        nosol = pd.DataFrame(nosolution, columns=["source Vdot", "Source Tin", "Source_Tout", "Sink Vdot", "Sink Tin",
                                                  "Sink Tout", "Power Set"])

        with pd.ExcelWriter(filename + ".xlsx") as writer:
            dfout.to_excel(writer, sheet_name="HP-Data")
            dfeff.to_excel(writer, sheet_name="COPH sorted")
            nosol.to_excel(writer, sheet_name="NoSolution")

        #   Resetting values
        self.sink_in.vdot = v0sink
        self.sink_out.vdot = self.sink_in.vdot
        self.source_in.vdot = v0source
        self.source_out.vdot = self.source_in.vdot
        return [dfout.head(10), dfeff.head(10), dfout.tail(10)]

    def operating_field_multi(self, mode, basefilename, source_temperatures, sink_temperatures):

        df_bestpower = pd.DataFrame([])
        df_lowpower = pd.DataFrame([])
        df_besteco = pd.DataFrame([])

        for source_temp in source_temperatures:
            for sink_temp in sink_temperatures:

                self.source_in.temperature = source_temp
                self.sink_out.temperature = sink_temp
                self.sink_in.temperature = sink_temp - 5.0

                filename = basefilename + str(source_temp) + "_" + str(sink_temp)
                out1 = self.operating_field_single(mode, sink_temp, filename)

                df_bestpower = df_bestpower.append(out1[0])
                df_besteco = df_besteco.append(out1[1])
                df_lowpower = df_lowpower.append(out1[2])

        with pd.ExcelWriter(basefilename + ".xlsx") as writer:
            df_bestpower.to_excel(writer, sheet_name="HP-BestPower")
            df_besteco.to_excel(writer, sheet_name="HP-BestEconomy")
            df_lowpower.to_excel(writer, sheet_name="HP-LowPower")

    def operating_field_fromFlow_single(self, target, filename):
        #  mode = 1 From_Flow    target = Vorlauftemperatur
        #  mode = 2 From Power   target = power in percent

        v0sink = self.sink_in.vdot
        v0source = self.source_in.vdot

        data = []
        nosolution = []

        for ivol_source in [-0.25, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5]:
            for ipower in [0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:

                self.source_in.vdot = v0source * (ivol_source + 1)
                self.source_out.vdot = self.source_in.vdot

                x = self.calculate_from_flow(ipower, target)

                if x != 0:
                    data.append([self.source_in.vdot, self.source_in.temperature, self.source_out.temperature,
                                 self.sink_in.vdot, self.sink_in.temperature, self.sink_out.temperature, ipower,
                                 self.TurboCore.power_condenser(),
                                 self.TurboCore.power_condenser() / self.TurboCore.power_electrical(),
                                 self.TurboCore.power_evaporator(), self.TurboCore.power_electrical(),
                                 self.TurboCore.temperature_suction(),
                                 self.TurboCore.temperature_discharge_saturated(),
                                 self.TurboCore.actual_values[0]
                                 ])
                    print("Power: ", self.TurboCore.power_condenser(),
                          self.TurboCore.power_condenser() / self.TurboCore.power_electrical())
                    self.write_log(x)
                else:
                    nosolution.append([self.source_in.vdot, self.source_in.temperature, self.source_out.temperature,
                                       self.sink_in.vdot, self.sink_in.temperature, self.sink_out.temperature,
                                       ipower])

        dfout = pd.DataFrame(data,
                             columns=["source Vdot", "Source Tin", "Source_Tout", "Sink Vdot", "Sink Tin", "Sink Tout",
                                      "Power Set", "Power Cond", "COPH", "Power Evap", "Power Electric", "Tsuction",
                                      " T Disch sat", "Actual Values"])
        dfout = dfout.sort_values(by='Power Cond', ascending=False)
        dfeff = dfout.sort_values(by='COPH', ascending=False)

        nosol = pd.DataFrame(nosolution, columns=["source Vdot", "Source Tin", "Source_Tout", "Sink Vdot", "Sink Tin",
                                                  "Sink Tout", "Power Set"])

        with pd.ExcelWriter(filename + ".xlsx") as writer:
            dfout.to_excel(writer, sheet_name="HP-Data")
            dfeff.to_excel(writer, sheet_name="COPH sorted")
            nosol.to_excel(writer, sheet_name="NoSolution")

        #   Resetting values
        self.sink_in.vdot = v0sink
        self.sink_out.vdot = self.sink_in.vdot
        self.source_in.vdot = v0source
        self.source_out.vdot = self.source_in.vdot
        return [dfout.head(10), dfeff.head(10), dfout.tail(10)]

    def operating_field_fromFLow_multi(self, basefilename, source_temperatures, sink_temperatures, deltaT_sink):

        df_bestpower = pd.DataFrame([])
        df_lowpower = pd.DataFrame([])
        df_besteco = pd.DataFrame([])

        for source_temp in source_temperatures:
            for sink_temp in sink_temperatures:
                self.source_in.temperature = source_temp
                self.sink_out.temperature = sink_temp
                self.sink_in.temperature = sink_temp - deltaT_sink

                filename = basefilename + str(source_temp) + "_" + str(sink_temp)
                out1 = self.operating_field_fromFlow_single(sink_temp, filename)

                df_bestpower = df_bestpower.append(out1[0])
                df_besteco = df_besteco.append(out1[1])
                df_lowpower = df_lowpower.append(out1[2])

        with pd.ExcelWriter(basefilename + ".xlsx") as writer:
            df_bestpower.to_excel(writer, sheet_name="HP-BestPower")
            df_besteco.to_excel(writer, sheet_name="HP-BestEconomy")
            df_lowpower.to_excel(writer, sheet_name="HP-LowPower")


class GasCoolerHeatPump(SimpleHeatPump):
    def __init__(self, compressor: TurboCor, condenser: Condenser, evaporator: Evaporator, source_in: Fluid,
                 source_out: Fluid, sink_in: Fluid, sink_out: Fluid, gascooler : GasCooler, hotwatersink_in : Fluid,
                 hotwaterwater_sink_out : Fluid):
        super().__init__(compressor, condenser, evaporator, source_in, source_out, sink_in,sink_out)
        # self.target_temperature = sink_out.temperature
        # self.TurboCore = compressor
        # self.condenser = condenser
        # self.evaporator = evaporator
        # self.source_in = source_in
        # self.source_out = source_out
        # self.sink_in = sink_in
        # self.sink_out = sink_out
        self.hotwater_sink_in = hotwatersink_in
        self.hotwater_sink_out = hotwaterwater_sink_out
        self.gascooler = gascooler
        self.p_gascooler = 0.0
        print("GasCoolerHeatPump initialized", compressor)

    def p_condenser(self):
        return self.TurboCore.power_condenser() - self.p_gascooler

    def calculate(self, frompower, t_suction, t_condensation):
        self.p_gascooler = self.gascooler.power(self.hotwater_sink_in, self.TurboCore.temperature_discharge())

        valid = self.TurboCore.calculate_frompower(frompower, t_suction, t_condensation)
        if valid != 1:
            print("Failure in Compressor valid = ", valid)
            quit("Failure")
        t_condensation_new = self.condenser.solvebalance(self.sink_in, self.p_condenser())[0]
        print("TurboCore Condenser Power: ", self.TurboCore.power_condenser(), t_condensation_new)
        self.sink_out.temperature = self.sink_in.heat(self.p_condenser() * 1000.0)
        print("Sink out Temperature: ", self.sink_out.temperature)
        t_suction_new = self.evaporator.solvebalance(self.source_in, -self.TurboCore.power_evaporator())[0]
        print("TurboCore Evaporator Power: ", self.TurboCore.power_evaporator(), t_suction_new)
        self.source_out.temperature = self.source_in.heat(-TurboCor.power_evaporator(self.TurboCore) * 1000.0)
        print("Temperature Source_out: ", self.source_out.temperature)
        return [frompower, t_suction_new, t_condensation_new, self.sink_out.temperature]

    def get_electrical_power(self):
        return self.TurboCore.power_electrical()



class AirWaterSimpleHeatpump(SimpleHeatPump):
    def __init__(self, compressor: TurboCor, condenser: Condenser, evaporator: Evaporator, source_in: Fluid,
                 source_out: Fluid, sink_in: Fluid, sink_out: Fluid, airunit: AirWaterunit):
        self.components = (sink_in, compressor, condenser, sink_out, source_in, airunit, evaporator, source_out)
        self.AirWaterunit = airunit
        super().__init__(compressor, condenser, evaporator, source_in, source_out, sink_in, sink_out)

    def get_electrical_power(self):
        return self.TurboCore.power_electrical() + self.AirWaterunit.Fan.electric_power

    def calculate(self, frompower, t_suction, t_condensation):   # Overriding calculate for AirWater HP

        q_ae = self.AirWaterunit.calculate_with_fan(self.AirWaterunit.Air, self.source_out)
        self.source_in.temperature = self.source_out.heat(q_ae)  # TODO super().calculate
        valid = self.TurboCore.calculate_frompower(frompower, t_suction, t_condensation)
        if valid != 1:
            print("Failure in Compressor valid = ", valid)
            quit("Failure")
        t_condensation_new = self.condenser.solvebalance(self.sink_in, self.p_condenser())[0]
        print("TurboCore Condenser Power: ", self.p_condenser(), t_condensation_new)
        self.sink_out.temperature = self.sink_in.heat(self.p_condenser() * 1000.0)
        print("Sink out Temperature: ", self.sink_out.temperature)
        t_suction_new = self.evaporator.solvebalance(self.source_in, -self.p_evaporator())[0]
        print("TurboCore Evaporator Power: ", self.p_evaporator(), t_suction_new)
        self.source_out.temperature = self.source_in.heat(-self.p_evaporator() * 1000.0)
        print("Temperature Source_out: ", self.source_out.temperature)
        return [frompower, t_suction_new, t_condensation_new, self.sink_out.temperature]

class AirWaterGasCoolerHeatpump(GasCoolerHeatPump):
    def __init__(self, compressor: TurboCor, condenser: Condenser, evaporator: Evaporator, source_in: Fluid,
                 source_out: Fluid, sink_in: Fluid, sink_out: Fluid, gascooler : GasCooler, hotwatersink_in: Fluid,
                 hotwaterwater_sink_out: Fluid, airunit: AirWaterunit):
        self.components = (sink_in, compressor, condenser, sink_out, source_in, airunit, evaporator, source_out)
        self.AirWaterunit = airunit
        super().__init__(compressor, condenser, evaporator, source_in, source_out, sink_in, sink_out, gascooler,
                         hotwatersink_in, hotwaterwater_sink_out)

    def get_electrical_power(self):
        return self.TurboCore.power_electrical() + self.AirWaterunit.Fan.electric_power

class TwoStageHeatPump(HeatPump):
    def __init__(self, hp1: HeatPump, hp2: HeatPump, compressor: Compressor, condenser: Condenser,
                 evaporator: Evaporator, source_in: Fluid, source_out: Fluid, sink_in: Fluid, sink_out: Fluid):
        super().__init__(compressor, condenser, evaporator, source_in, source_out, sink_in, sink_out)
        self.components = (hp1.components + hp2.components)
