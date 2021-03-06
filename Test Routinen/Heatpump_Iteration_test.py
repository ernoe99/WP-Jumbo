import os as os
import time

from component_use import *

# Initialize Heatpump

def output_organize(base):
    command = "mkdir output\\" + base
    rv = os.system(command)
    command = "move " + base + "*.xlsx " + "output\\" + base
    rv = os.system(command)
    command = "move *.log " + "output\\" + base
    rv = os.system(command)
    command = "move *.txt " + "output\\" + base
    rv = os.system(command)


# Variables


x = Thermalia_375.calculate_from_power(50.0)

print("RESULT: ", x)

y = Thermalia375_Basic.calculate_from_power(50.0)

print("RESULT: ", y)

ye = Thermalia375_Basic_Eco.calculate_from_power(50.0)

print("RESULT: ", ye)

HW55_RL.temperature = 25.0

z = Thermalia375_Basic_Eco_W.calculate_from_power(50.0)

print("RESULT: ", z)

Thermalia375_Basic_Eco_W.hotwater_sink_in.vdot = 0.0
Thermalia375_Basic_Eco_W.hotwater_sink_out.vdot = 0.0

zw = Thermalia375_Basic_Eco_W.calculate_from_power(50.0)

print("RESULT: ", zw)

zw = Thermalia375_Basic_Eco_W0.calculate_from_power(50.0)

print("RESULT: ", zw)


# Thermalia_375F.source_in.vdot = 12.5
# Thermalia_375F.sink_in.vdot = 13.65212997
# Thermalia_375.source_out.vdot = Thermalia_285.source_in.vdot
# Thermalia_375.sink_out.vdot = Thermalia_375.sink_in.vdot
# Thermalia_375F.source_in.temperature = - 10.0
#
# x = Thermalia_375.calculate_from_flow(90, 35.0)

# source_temperature = -5.0
#
# Thermalia_375.source_in.temperature = source_temperature
# Thermalia_375F.source_in.temperature = source_temperature
#
#
# df = pd.DataFrame([])
# cf = pd.DataFrame([])
#
#
# out1 = Thermalia_375.operating_field_single(1, 35.0, "Thermalia_375_fromFlow35")
# # out2 = Thermalia_375F.operating_field_single(1, 35.0, "Thermalia_375F_fromFlow35")
#
# df = df.append(out1[0])
# print(df)
# cf = cf.append(out1[1])
# print(cf)
# cf = cf.append(out1[0])
# print(cf)
#
sctemp = []
for item in range(-14, 21, 1):
    sctemp.append(float(item))
sitemp = []
for item in range(10, 76, 1):
    sitemp.append(float(item))

print(sctemp)
print(sitemp)

dtsink = 5.0

start_time = time.time()

# base = "CCC_Thermalia375_Basic"
# Thermalia375_Basic.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
#                                                             sink_temperatures=sitemp, deltaT_sink=dtsink)
# output_organize(base)
#
# base = "CCC_Thermalia375_Basic_Eco"
# Thermalia375_Basic_Eco.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
#                                                            sink_temperatures=sitemp, deltaT_sink=dtsink)
# output_organize(base)
#
# base = "CCC_Thermalia375_Basic_Eco_W"
# Thermalia375_Basic_Eco_W.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
#                                                             sink_temperatures=sitemp, deltaT_sink=dtsink)
# output_organize(base)
#
# base = "CCC_Thermalia375_Performance"
# Thermalia375_Performance.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
#                                                             sink_temperatures=sitemp, deltaT_sink=dtsink)
# output_organize(base)
#
#
# base = "CCC_Thermalia375_Basic_Eco_W0"
# Thermalia375_Basic_Eco_W0.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
#                                                             sink_temperatures=sitemp, deltaT_sink=dtsink)
# output_organize(base)


######################################### TGH 285 ################################################
base = "CCC_Thermalia285_Basic"
Thermalia285_Basic.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
                                                            sink_temperatures=sitemp, deltaT_sink=dtsink)
output_organize(base)

base = "CCC_Thermalia285_Basic_Eco"
Thermalia285_Basic_Eco.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
                                                           sink_temperatures=sitemp, deltaT_sink=dtsink)
output_organize(base)

base = "CCC_Thermalia285_Basic_Eco_W"
Thermalia285_Basic_Eco_W.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
                                                            sink_temperatures=sitemp, deltaT_sink=dtsink)
output_organize(base)

base = "CCC_Thermalia285_Performance"
Thermalia285_Performance.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
                                                            sink_temperatures=sitemp, deltaT_sink=dtsink)
output_organize(base)


base = "CCC_Thermalia285_Basic_Eco_W0"
Thermalia285_Basic_Eco_W0.operating_field_fromFLow_multi(basefilename=base, source_temperatures=sctemp,
                                                            sink_temperatures=sitemp, deltaT_sink=dtsink)
output_organize(base)

print("--- %s seconds ---" % (time.time() - start_time))

#

# for ipower in range(10, 11):
#
#     with open('logfile.txt', 'a') as log:
#         log.write("Power: " + str(ipower * 10.0) + "\n")
#
#     val_iter = Thermalia_285.calculate(t_target=requested_T, frompower=ipower * 10.0, t_suction=t_suc,
#                                        t_condensation=t_cond)
#     print(" New Values", val_iter[0], val_iter[1], val_iter[2] + (val_iter[3] - val_iter[4]))
#     print("Iteration Values: ", val_iter)
#     with open('logfile.txt', 'a') as log:
#         log.write(
#             " New Values {0}  {1}  {2}\n".format(str(val_iter[0]), str(val_iter[1]),
#                                                  str(val_iter[2] + (val_iter[3] - val_iter[4]))))
#         log.write("Iteration Values: " + str(val_iter) + "\n")
#
#     for i in range(1, 50):
#         val_iter = Thermalia_285.calculate(t_target=requested_T, frompower=ipower * 10.0, t_suction=val_iter[1],
#                                            t_condensation=val_iter[2] + (val_iter[3] - val_iter[4]))
#         print(" New Values", val_iter[0], val_iter[1], val_iter[2] + (val_iter[3] - val_iter[4]))
#         print("Iteration Values: ", val_iter)
#         with open('logfile.txt', 'a') as log:
#             log.write(
#                 " New Values {0}  {1}  {2}\n".format(str(val_iter[0]), str(val_iter[1]),
#                                                      str(val_iter[2] + (val_iter[3] - val_iter[4]))))
#             log.write("Iteration Values: " + str(val_iter) + "\n")
#             log.write("Power Data from Compressor:{0} {1} {2} {3} {4}\n".format(str(TTH375.power_electrical()),
#                                                                                 str(TTH375.power_evaporator()),
#                                                                                 str(TTH375.power_condenser()),
#                                                                                 str(TTH375.power_economizer()),
#                                                                                 str(TTH375.power_condenser() /
#                                                                                     TTH375.power_electrical())))
#
# print(TTH375.actual_values)
