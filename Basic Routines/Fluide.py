import CoolProp.CoolProp as CP


class Refrigerant:
    def __init__(self, name='R32', temp=20.0, mdot=1.0):
        self.name = name
        t1 = 273.15 + temp
        p1 = CP.PropsSI('P', 'T', t1, 'Q', 0.5, self.name)
        h1x = CP.PropsSI('H', 'T', t1, 'Q', 0.5, self.name)

        self.state = [mdot, t1, p1, h1x]


r32 = Refrigerant('R32')
print(r32.state)


fluid = 'R32'

pressure_at_critical_point = CP.PropsSI(fluid,'pcrit')

print(pressure_at_critical_point/1.0E5)
#
h1 = CP.PropsSI('H', 'T', 273.15, 'P', 5.0E5, fluid)

print(h1)
H_V = CP.PropsSI('H', 'P', 101325, 'Q', 1, fluid);
print(H_V)
