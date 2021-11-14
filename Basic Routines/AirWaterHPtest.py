# Air Water Heat pump test driver

import Fluids
import HeatExchanger
import PumpsandFans
import Heatpump
from component_use import W35_RL_350, W35_350

outdoor_air = Fluids.WetAir(2.0, 10.0, 0.9)

ae_375fan = PumpsandFans.Fan(0.7, 20.0, 70.0)

ae_375 = HeatExchanger.AirWaterunit(ae_375fan, outdoor_air, 5000.0, 100.0)

s_in = Fluids.Brine(-1.0, 17.6)

x = ae_375.calculate(outdoor_air, s_in, 2)

print(x)

ae_375.set_air_volume(20.0, s_in)

print(ae_375.Air.vdot)
x = ae_375.calculate(outdoor_air, s_in, 2)

print(x)

s_out = s_in
s_out.temperature = s_in.heat(x)

print(s_out.temperature)

Bellaria375_Performance = Heatpump.AirWaterSimpleHeatpump(compressor=TTH375, condenser=B427Hx400, evaporator=CSF2880,
                                                   sink_in=W35_RL_350, sink_out=W35_350, source_in=s_in,
                                                   source_out=s_out, airunit=ae_375)







