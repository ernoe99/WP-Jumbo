# This is the file for usable components - bottom up
import pandas as pd

import Fluids
import HeatExchanger
import Heatpump
import PumpsandFans

from TurboCor import TTH375, TGH285, TGH285_noEcon, TTH375_noEcon

# Heat Exchangers

# Evaporators

# PROVIDES

#  CSF2880_fitted = Evaporator_fitted(306, 0.4183, 250, (-13.3, 24.01, -0.2405))  altes fitting
#  CSF2880_fitted = Evaporator_fitted(306, 0.4183, 250, (-22.616, 110.69, -3.007))
CSF2880 = HeatExchanger.Evaporator_fitted(306, 0.4183, 250, (404.92, 45.341, -0.6457))


# Condensers

# SWEP


B427Hx140 = HeatExchanger.Condenser_fitted(140, 0.192, 750.0, (1058, 64.855, 0.0))
B427Hx400 = HeatExchanger.Condenser_fitted(400, 0.194, 750.0, (719.56, 37.262, 0))

B65Hx348 = HeatExchanger.Evaporator_fitted(348, 0.257, 900, (100.75, 78.439, -1.6407))


# Pumps

PumpSource = PumpsandFans.Pump(20.0, 2.0, 0.8)   # Dummy Source pump
PumpSink = PumpsandFans.Pump(20.0, 2.0, 0.8)   # Dummy Sink pump
PumpHW = PumpsandFans.Pump(10.0, 1.0, 0.8)   # Dummy heatwater pump


# Heaters


# Air water units

# Gas cooler

GCdummy = HeatExchanger.GasCooler(PumpHW, 1000.0)



# Standard Sources

B0_Geo = Fluids.Brine(0, 25.0)
source_out = Fluids.Brine(-3, 25.0)
W35_250 = Fluids.Heatwater(35.0, 11.96)
W35_350 = Fluids.Heatwater(35.0, 16.75)
W35_RL_250 = Fluids.Heatwater(30.0, 11.96)
W35_RL_350 = Fluids.Heatwater(30.0, 16.75)
HW55_RL = Fluids.Heatwater(55.0, 3.0)
HW57_VL = Fluids.Heatwater(57.0, 3.0)
HW55_RL_0 = Fluids.Heatwater(55.0, 0.0)   # Fuer Heisswasser abgeschaltet
HW57_VL_0 = Fluids.Heatwater(57.0, 0.0)   # Fuer Heisswasser abgeschaltet

# Heat Pumps
Thermalia_375 = Heatpump.SimpleHeatPump(compressor=TTH375, condenser=B427Hx140, evaporator=CSF2880,
                                        sink_in=W35_RL_350, sink_out=W35_350, source_in=B0_Geo, source_out=source_out)
Thermalia_375F = Heatpump.SimpleHeatPump(compressor=TTH375, condenser=B427Hx400, evaporator=CSF2880,
                                         sink_in=W35_RL_350, sink_out=W35_350, source_in=B0_Geo, source_out=source_out)
Thermalia_375E = Heatpump.SimpleHeatPump(compressor=TTH375, condenser=B427Hx140, evaporator=B65Hx348,
                                         sink_in=W35_RL_350, sink_out=W35_350, source_in=B0_Geo, source_out=source_out)


Thermalia_285 = Heatpump.SimpleHeatPump(compressor=TGH285, condenser=B427Hx140, evaporator=CSF2880,
                                        sink_in=W35_RL_250, sink_out=W35_250, source_in=B0_Geo, source_out=source_out)
Thermalia_285F = Heatpump.SimpleHeatPump(compressor=TGH285, condenser=B427Hx400, evaporator=CSF2880,
                                         sink_in=W35_RL_250, sink_out=W35_250, source_in=B0_Geo, source_out=source_out)

print("Heat pumps first stage initiated")

# Heat pumps - second stage

Thermalia375_Basic = Heatpump.SimpleHeatPump(compressor=TTH375_noEcon, condenser=B427Hx140, evaporator=B65Hx348,
                                             sink_in=W35_RL_350, sink_out=W35_350, source_in=B0_Geo,
                                             source_out=source_out)

Thermalia375_Basic_Eco = Heatpump.SimpleHeatPump(compressor=TTH375, condenser=B427Hx140, evaporator=B65Hx348,
                                                 sink_in=W35_RL_350, sink_out=W35_350, source_in=B0_Geo,
                                                 source_out=source_out)

Thermalia375_Basic_Eco_W = Heatpump.GasCoolerHeatPump(compressor=TTH375, condenser=B427Hx140, evaporator=B65Hx348,
                                                      sink_in=W35_RL_350, sink_out=W35_350, source_in=B0_Geo,
                                                      source_out=source_out,
                                                      gascooler=GCdummy, hotwatersink_in=HW55_RL,
                                                      hotwaterwater_sink_out=HW57_VL)
# Heatpump with gascooler but no water flow
Thermalia375_Basic_Eco_W0 = Heatpump.GasCoolerHeatPump(compressor=TTH375, condenser=B427Hx140, evaporator=B65Hx348,
                                                       sink_in=W35_RL_350, sink_out=W35_350, source_in=B0_Geo,
                                                       source_out=source_out,
                                                       gascooler=GCdummy, hotwatersink_in=HW55_RL_0,
                                                       hotwaterwater_sink_out=HW57_VL_0)

Thermalia375_Performance = Heatpump.SimpleHeatPump(compressor=TTH375, condenser=B427Hx400, evaporator=CSF2880,
                                                   sink_in=W35_RL_350, sink_out=W35_350, source_in=B0_Geo,
                                                   source_out=source_out)

############################################################  TTH 285 ##########################
Thermalia285_Basic = Heatpump.SimpleHeatPump(compressor=TGH285_noEcon, condenser=B427Hx140, evaporator=B65Hx348,
                                             sink_in=W35_RL_250, sink_out=W35_250, source_in=B0_Geo,
                                             source_out=source_out)

Thermalia285_Basic_Eco = Heatpump.SimpleHeatPump(compressor=TGH285, condenser=B427Hx140, evaporator=B65Hx348,
                                                 sink_in=W35_RL_250, sink_out=W35_250, source_in=B0_Geo,
                                                 source_out=source_out)

Thermalia285_Basic_Eco_W = Heatpump.GasCoolerHeatPump(compressor=TGH285, condenser=B427Hx140, evaporator=B65Hx348,
                                                      sink_in=W35_RL_250, sink_out=W35_250, source_in=B0_Geo,
                                                      source_out=source_out,
                                                      gascooler=GCdummy, hotwatersink_in=HW55_RL,
                                                      hotwaterwater_sink_out=HW57_VL)
# Heatpump with gascooler but no water flow
Thermalia285_Basic_Eco_W0 = Heatpump.GasCoolerHeatPump(compressor=TGH285, condenser=B427Hx140, evaporator=B65Hx348,
                                                       sink_in=W35_RL_250, sink_out=W35_250, source_in=B0_Geo,
                                                       source_out=source_out,
                                                       gascooler=GCdummy, hotwatersink_in=HW55_RL_0,
                                                       hotwaterwater_sink_out=HW57_VL_0)

Thermalia285_Performance = Heatpump.SimpleHeatPump(compressor=TGH285, condenser=B427Hx400, evaporator=CSF2880,
                                                   sink_in=W35_RL_250, sink_out=W35_250, source_in=B0_Geo,
                                                   source_out=source_out)


#######################################

#  STARTING TESTS

#######################################

def test_loop_hp(hp, t_source_in, vdot_source_in, t_rl, vdot_sink_in, filename):


    # hp = Thermalia_285

    x = hp.calculate_from_power(80.0)

    v0sink = vdot_sink_in
    v0source = vdot_source_in

    hp.source_in.vdot = vdot_source_in
    hp.source_out.vdot = hp.source_in.vdot
    hp.source_in.temperature = t_source_in
    hp.sink_in.vdot = vdot_sink_in
    hp.sink_out.vdot = hp.sink_in.vdot
    hp.sink_in.temperature = t_rl
    hp.sink_out.temperature = t_rl + 5  # nicht gebraucht - nur Info

    #   x = hp.calculate_from_power(ipower)

    data = []
    nosolution = []

    for ivol_sink in [-0.5, -0.25, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5]:
        for ivol_source in [-0.5, -0.25, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5]:
            for ipower in [50, 60, 70, 80, 90, 100]:

                hp.sink_in.vdot = v0sink * (ivol_sink + 1)
                hp.sink_out.vdot = hp.sink_in.vdot
                hp.source_in.vdot = v0source * (ivol_source + 1)
                hp.source_out.vdot = hp.source_in.vdot

                x = hp.calculate_from_power(ipower)

                print(v0sink * (ivol_sink + 1), v0source * (ivol_source + 1), ipower)

                if x != 0:
                    data.append([hp.source_in.vdot, hp.source_in.temperature, hp.source_out.temperature,
                                 hp.sink_in.vdot, hp.sink_in.temperature, hp.sink_out.temperature, ipower,
                                 hp.TurboCore.power_condenser(), hp.TurboCore.power_condenser() / hp.TurboCore.power_electrical(),
                                 hp.TurboCore.power_evaporator(), hp.TurboCore.power_electrical(),
                                 hp.TurboCore.temperature_suction(), hp.TurboCore.temperature_discharge_saturated(),
                                 hp.TurboCore.actual_values[0]
                                 ])
                    print("Power: ", hp.TurboCore.power_condenser(), hp.TurboCore.power_condenser() / hp.TurboCore.power_electrical())
                    hp.write_log(x)
                else:
                    nosolution.append([hp.source_in.vdot, hp.source_in.temperature, hp.source_out.temperature,
                                 hp.sink_in.vdot, hp.sink_in.temperature, hp.sink_out.temperature, ipower])

    dfout = pd.DataFrame(data, columns=["source Vdot", "Source Tin", "Source_Tout", "Sink Vdot", "Sink Tin", "Sink Tout",
                                        "Power Set", "Power Cond", " COPH", "Power Evap", "Power Electric", "Tsuction",
                                        " T Disch sat", "Actual Values"])

    nosol = pd.DataFrame(nosolution, columns=["source Vdot", "Source Tin", "Source_Tout", "Sink Vdot", "Sink Tin",
                                              "Sink Tout", "Power Set"])

    with pd.ExcelWriter(filename + ".xlsx") as writer:
        dfout.to_excel(writer, sheet_name="HP-Data")
        nosol.to_excel(writer, sheet_name="NoSolution")


# start_time = time.time()
#
# test_loop_hp(Thermalia_285, 0.0, 20.0, 25.0, 11.96, 'Thermalia__285_A')
#
# test_loop_hp(Thermalia_285F,  0.0, 20.0, 25.0, 11.96, 'Thermalia__285F_A')
#
# test_loop_hp(Thermalia_375, 0.0, 30.0, 25.0, 18, 'Thermalia__375_A')
#
# test_loop_hp(Thermalia_375F, 0.0, 30.0, 25.0, 18, 'Thermalia__375_A')
#
# print("--- %s seconds ---" % (time.time() - start_time))
