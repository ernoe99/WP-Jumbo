import math


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
        return [self.mdot(self.temperature), self.vdot, self.temperature, self.density(self.temperature),
                self.enthalpy(), self.cp]

    def capacity(self, temp):
        return self.mdot(temp) * self.cp

    def heat(self, qdot):
        return self.temperature + qdot / self.capacity(self.temperature)

class Heatwater(Fluid):
    def __init__(self, temp, volflow):
        super().__init__(temp, volflow)
        self.cp = 4184.0  # override cp for Heatwater

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
        """

        :rtype: object
        """
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

        return K3 * ((K2 * temp) / (K3 + temp) + math.log(phi)) / (K2 * K3 / (K3 + temp) - math.log(phi))

    def t_h_x_p(self, h, x, p):   # Umkehrfunktion zu h_t_x_p
        return (h - x * self.hv_t(0.0)) / (1.006 + 1.86*x)

    def phi_t_x_p(self, temp, xact, p):
        pds = self.pds_t(temp)
        pd = p * xact / (0.622 + xact)
        return min(pd / pds, 1.0)

    def h_t_x_p(self, temp, xact, p):
        xs = self.xs_t_p(temp, p)
        cpx = 1.006
        if xact > xs:
            h = 1000.0*(temp * cpx + xs * (self.hv_t(0) + 1.86 * temp))
            xliquid = xact - xs
        else:
            h = 1000*(temp * cpx + xact * (self.hv_t(0) + 1.86 * temp))
            xliquid = 0
        return [h, xliquid]

    def condens(self, t_target, x, p):
        xs = self.xs_t_p(t_target, p)
        return max((x - xs), 0.0)


class Burned_Air:
    def __init__(self, temperature, mdot):
        self.temperature = temperature
        self.mdot = mdot
        self.cp = 1006.0

    def capacity(self, temp):
        return self.mdot * self.cp

    def heat(self, qdot):
        return self.temperature + qdot / self.capacity(1)

class HAir(Air):
    def __init__(self, temp, volflow, rf):
        super().__init__(temp, volflow)
        self.cp = 1006.0
        self.rF = rf     # relative Feuchte o-1
        self.pressure = 101325.0
        self.Ktemperature = temp + 273.15
        self.temperature = temp




# Initialisierung wichtiger stati

air_2 = WetAir(2.0, 1000.0, 0.8)  # Achtung luft in m3/s Wasser in l/s
ground_0 = Brine(0.0, 25.0)       # 0 Grad 25 l/s Wasser Glykol 60/40

h = air_2.h_t_x_p(2, air_2.x_t_phi_p(2, 0.8, 101325), 101325)
print(h)

h = air_2.h_t_x_p(-3, air_2.x_t_phi_p(-3, 1, 101325), 101325)
print(h)

print(air_2.hv_t(0))

for i in [3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7]:
    t = float(i)
    xs = air_2.xs_t_p(t,101325)
    hx = air_2.h_t_x_p(t,xs,101325)
    # print(f"{t}  {xs} {hx[0]}")
    print(f"{hx[0]} ")
