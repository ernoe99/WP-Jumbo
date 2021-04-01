from Fluids import Brine, Heatwater, Fluid
from HeatExchanger import Evaporator_fitted, Evaporator, Condenser_fitted, Condenser, HeatExchanger, TwophaseHX
import pandas as pd


def test_condenser(condenser, t_sink_in, vdot_sink_in, pcondenser):
    fluid_in = Heatwater(t_sink_in, vdot_sink_in)

    return condenser.solvebalance(fluid_in, pcondenser)


def loop_test_condenser(condenser,t_sink_in, vdot_sink_in, pcondenser):
    data = []
    for iflow in [-0.5, -0.25, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.25, 0.5]:

        vdotx = vdot_sink_in * (1 + iflow)

        data.append([t_sink_in, vdotx, pcondenser, test_condenser(condenser, t_sink_in, vdotx, pcondenser)])

        print(data)


def test_evaporator(evaporator, t_source_in, vdot_source_in, pevaporator):
    fluid_in = Heatwater(t_source_in, vdot_source_in)

    return evaporator.solvebalance(fluid_in, pevaporator)


def loop_test_evaporator(evaporator, t_source_in, vdot_source_in, pevaporator):
    data = []
    for iflow in [-0.5, -0.25, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.25, 0.5]:

        vdotx = vdot_source_in * (1 + iflow)

        data.append([t_source_in, vdotx, pevaporator, test_evaporator(evaporator, t_source_in, vdotx, pevaporator)[0],
                                                                      t_source_in + pevaporator / 4.18 / vdotx])

    df = pd.DataFrame(data, columns=["Tsource", "vdot source", " Power", "T2phase", " Toutsource"])
    print(df)



#####  ALFA LAVAL  #####
AL_CB210_134AH_F = Condenser(134, 0.2485, 972.0)  # fitted

# Fitting of Alfa Laval Kondensator
# xf = Heatwater(28, 51652.0/3600)
# ufit = AL_KB210_134AH_F.fithx2phase(xf, 300.0, 35.0)
# print("UFIT: ", ufit)
# print("NEW Temperature: ", xf.heat(300000.0))

hf = Brine(0.0, 61134 / 3600.0)
AL_ACH1000DQ_306AH = Evaporator(306, 0.4183, 233.0)  # fitted
hfit = AL_ACH1000DQ_306AH.fithx2phase(hf, -250.0, -6.5)

print("HFIT: ", hfit)
print("NEW Temperature: ", hf.heat(-250000.0))

##### P R O V I D E S ########

vol11 = 275000.0 / 3600.0 / (12 - 7)
print("vol11: ", vol11)

cpx = 4180
f11 = Heatwater(12, vol11)
CSF2880 = Evaporator(306, 0.4183, 250)


f11.cp = cpx
p11 = -180.0
f11.vdot = -p11 * 1000.0 / cpx / 5.0

hfit = CSF2880.fithx2phase(f11, p11, 6.56)  # Provides Evaporator
print("HFIT: ", hfit)
print("NEW Temperature: ", f11.heat(p11 * 1000.0), f11.vdot)

p11 = -220.0
f11.vdot = -p11 * 1000.0 / cpx / 5.0

hfit = CSF2880.fithx2phase(f11, p11, 6.48)  # Provides Evaporator
print("HFIT: ", hfit)
print("NEW Temperature: ", f11.heat(p11 * 1000.0), f11.vdot)

p11 = -250.0
f11.vdot = -p11 * 1000.0 / cpx / 5.0

hfit = CSF2880.fithx2phase(f11, p11, 6.40)  # Provides Evaporator
print("HFIT: ", hfit)
print("NEW Temperature: ", f11.heat(p11 * 1000.0), f11.vdot)

p11 = -275.0
f11.vdot = -p11 * 1000.0 / cpx / 5.0
hfit = CSF2880.fithx2phase(f11, p11, 6.32)  # Provides Evaporator
print("HFIT: ", hfit)
print("NEW Temperature: ", f11.heat(p11 * 1000.0), f11.vdot)

p11 = -320.0
f11.vdot = -p11 * 1000.0 / cpx / 5.0
hfit = CSF2880.fithx2phase(f11, p11, 6.16)  # Provides Evaporator
print("HFIT: ", hfit)
print("NEW Temperature: ", f11.heat(p11 * 1000.0), f11.vdot)

f12 = Fluid(26.0, 17.5)
CCF2880315_6P = Condenser(134, 0.2485, 937.0)  # Provides Condenser
print("CCF2880315_6P.ua: ", CCF2880315_6P.ua)

CSF2880_fitted = Evaporator_fitted(306, 0.4183, 250, (-22.616, 110.69, -3.007))

print("stop")
print(" Fitted: ", CSF2880_fitted.ua, CSF2880_fitted.ua, CSF2880_fitted.ua)

print(" Fitted: ", CSF2880_fitted.uaf(9.0), CSF2880_fitted.uaf(11.0), CSF2880_fitted.uaf(21.0))
print(" Fitted: ", CSF2880_fitted.ntu(f11))

loop_test_evaporator(CSF2880_fitted, 0, 11.96, -250.0)

f13 = Heatwater(0, 11.96)
CSF2880_fitted.loop_vdot_HX(f13, -250.0)

#SWEP


B427Hx140 = Condenser(140, 0.192, 1500.0)

power = [150, 200, 250, 300, 350, 400.0]
volu = [7.18, 9.573, 11.97, 14.36, 16.75, 19.14]

for i in range(0, len(volu)):
    f12 = Heatwater(temp=30, volflow=volu[i])
    hfit = B427Hx140.fithx2phase(f12, power[i], 37.0)
    print(" Data: ", volu[i], power[i], hfit)

B427Hx140 = Condenser_fitted(140, 0.192, 750.0, (1058, 64.855, 0.0))

area = 26.7
print(" Fitted: ", B427Hx140 .uaf(9.572)/area, B427Hx140 .uaf(11.97)/area, B427Hx140 .uaf(19.14)/area)

f13 = Heatwater(temp=38, volflow=14.36)
B427Hx140.loop_vdot_HX(f13, 300)

print(B427Hx140.solvebalance(f13, 300.0))




B427Hx300 = Condenser(300, 0.194, 596.0)

# power = [150, 200, 250, 300, 350, 400]
# volu = [7.18, 9.573, 11.97, 14.36, 16.75, 19.15]
#
# for i in range(0, 6):
#     f12 = Fluid(temp=30, volflow=volu[i])
#     hfit = B427Hx300.fithx2phase(f12, power[i], 37.0)
#     print(" Data: ", volu[i], power[i], hfit)
#
# B427Hx300_fitted = Condenser_fitted(300, 0.194, 750.0, (-0.3068, 38.37,  -0.0024))
#
# print(" Fitted: ", B427Hx300_fitted.uaf(9.0), B427Hx300_fitted.uaf(11.0), B427Hx300_fitted.uaf(21.0))
#
# for i in range(0, 6):
#     f12 = Fluid(temp=50, volflow=volu[i])
#     hfit = B427Hx300.fithx2phase(f12, power[i], 57.0)
#     print(" Data: ", volu[i], power[i], hfit)

B427Hx400 = Condenser_fitted(400, 0.194, 750.0, (719.56, 37.262, 0))


f13 = Heatwater(temp=38, volflow=14.36)
B427Hx400.loop_vdot_HX(f13, 300)

print(B427Hx400.solvebalance(f13, 300.0))

# B427Hx120 = Condenser_fitted(100, 0.194, 750.0, (37.262/4.0, 719.56/4.0, 0))
# B427Hx120.loop_vdot_HX(f13, 350)