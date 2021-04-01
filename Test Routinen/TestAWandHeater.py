from Fluids import *
from HeatExchanger import *

Ultragas500 = Heater(500.0, 0.98, 0.09, 1.0, 750.0)

hotflow = Ultragas500.hot_flow(500.0, 1.0)

water = Heatwater(35.0, 12.0)
print(water.heat(Ultragas500.power(1.0)))
power_water = Ultragas500.power_water(hotflow, water)
print(Ultragas500.power_water(hotflow, water))

rest_heat = Ultragas500.power(1) - power_water

print("***", power_water, rest_heat)

restgas = Burned_Air(0, hotflow.mdot)
restgas.temperature = restgas.heat(rest_heat)
print(restgas.temperature)

x = Ultragas500.exhaust_gas_temperature(1.0)
print(Ultragas500.exhaust_gas_temperature(1.0))

fan1 = Fan(0.72, 100.0, 200.0)
print(fan1.electric_power)

fraction = 0.8  # actor
air = WetAir(2.0, fan1.max_volume * fraction, 0.5)  # reactor
water = Brine(-2.0, 25.0)

print("FAn1.calculate ", fan1.calculate(air))
print("Fan1.electric Power: ", fan1.electric_power)

ae1 = AirWaterunit(fan1, 10.0, 40000.0)

x = ae1.calculate(air, water, 2)
print(x)
x = ae1.calculate(air, water, 1)
print(x)

new_t = water.heat(x)
print(new_t)

t = ae1.calculate_with_fan(air, water)
print(t)

print("Electric Power Fan: ", ae1.Fan.electric_power, air.vdot, ae1.Fan.d_pressure(air))



