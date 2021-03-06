import numpy as np
import pandas as pd
from scipy.interpolate import griddata


from Fluids import Fluid
from PumpsandFans import Pump


class Compressor:
    def __init__(self, hsuc, hdisch, hliquid, mflow):
        self.enthalpy_liquid_subcooled = hliquid
        self.enthalpy_suction = hsuc
        self.enthalpy_discharge = hdisch
        self.massflow_discharge = mflow

    def q_condenser(self):
        return (self.enthalpy_discharge - self.enthalpy_liquid_subcooled) * self.massflow_discharge


class TurboCor(Compressor):
    def __init__(self, fname: str):
        # optimal data w/o Economizer
        fname = "..\\Datenfiles\\" + fname
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

        valid = self.calculate_frompower(frompower=50.0, t_suction=0.0, t_condensation=37.5)
        print("TurboCor Init: ", fname, "Minimum Power 50% 0 37.5 Grad: ", self.power_minimum(),
              " Maximum Power: ", self.power_maximum())
        self.economizer = 0

    def calculate_0power_econ(self, t_suction, t_condensation):

        idx = int((t_suction + 18) / 3)
        jdx = int((t_condensation + 5) / 2.5)
        vec0 = np.empty(38)
        valid = 0

        if -17.99 <= t_suction <= 30 and -5 <= t_condensation <= 69.99:
            rs = self.zeropower_econ_data

            if rs[idx, jdx, 14] != 0 and rs[idx + 1, jdx, 14] != 0 and rs[idx, jdx + 1, 14] != 0 and rs[
                idx + 1, jdx + 1, 14] != 0:

                for i in range(1, 39, 1):  # linear interpolation
                    w11 = rs[idx, jdx, i] + (rs[idx + 1, jdx, i] - rs[idx, jdx, i]) / 3.0 * (
                                t_suction - rs[idx, jdx, 1])
                    w12 = rs[idx, jdx + 1, i] + (rs[idx + 1, jdx + 1, i] - rs[idx, jdx + 1, i]) / 3.0 * (
                            t_suction - rs[idx, jdx + 1, 1])
                    w23 = w11 + (w12 - w11) / 2.5 * (t_condensation - rs[idx, jdx, 2])
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
        valid = self.calculate_0power_econ(t_suction, t_condensation)
        if valid == 0:
            valid = self.calculate_0power_noecon(t_suction, t_condensation)
        return valid  # Results in self.actual_values_0Power

    def get_values(self, t_suction, t_condensation, rs):
        id = int((t_suction + 18) / 3)
        jd = int((t_condensation + 5) / 2.5)
        vec0 = np.empty(38)

        if -17.99 <= t_suction <= 30 and -5 <= t_condensation <= 69.99:

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
        return self.actual_values[0] + self.actual_values[20]  # +[20] ist superheat

    def temperature_discharge_saturated(self):
        return self.actual_values[19]

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

    def power_minimum(self):
        return self.actual_values[36]

    def power_maximum(self):
        return self.actual_values[37]


class TurboCor_noEcon(TurboCor):
    def __init__(self, fname: str):
        super(TurboCor_noEcon, self).__init__(fname)

    def calculate_0power(self, t_suction, t_condensation):
        valid = self.calculate_0power_noecon(t_suction, t_condensation)
        return valid  # Results in self.actual_values_0Power

    def calculate_zeropower(self, t_suction, t_condensation):

        valid = self.get_values(t_suction, t_condensation, self.zeropower_noecon_data)
        self.economizer = 0
        return valid

    def calculate_frompower(self, frompower, t_suction, t_condensation):

        id1 = int(frompower / 25.0)
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


class PolyScroll(Compressor):
    def __init__(self, fname: str, rpm):
        self.rpm = rpm

        # Compressor data from Polynoms with speeds
        fname = "..\\Datenfiles\\" + fname
        fcsv = fname
        pddata = pd.read_csv(fcsv, sep=';')
        linedat = pddata.to_numpy(copy=True)
        self.poly_data = linedat

    def calculate_fromrpm(self, speed, t_suction, t_condensation):
        data = self.poly_data
        ishift = 7
        rpm = self.rpm
        capacity = np.empty(4, float)  # nur fuer die Dimension der Liste
        power = np.empty(4, float)    # nur fuer die Dimension der Liste
        flow_rate =np.empty(4, float)     # nur fuer die Dimension der Liste
        current = np.empty(4, float)   # nur fuer die Dimension der Liste

        for isp in range(0, len(self.rpm)):
            ic = ishift * isp + 1
            capacity[isp] = (data[0, ic] + t_suction * data[1, ic] + t_condensation * data[2, ic] +
                            t_suction**2 * data[3, ic] + t_suction * t_condensation * data[4, ic] +
                            t_condensation**2 * data[5, ic] + t_suction**3 *
                            data[6, ic] + t_suction**2 * t_condensation * data[7, ic] +
                            t_suction * t_condensation**2 * data[8, ic] +
                            t_condensation**3 * data[9, ic] )
            ic += 1
            power[isp] = (data[0, ic] + t_suction * data[1, ic] + t_condensation * data[2, ic] +
                            t_suction**2 * data[3, ic] + t_suction * t_condensation * data[4, ic] +
                            t_condensation**2 * data[5, ic] + t_suction**3 *
                            data[6, ic] + t_suction**2 * t_condensation * data[7, ic] +
                            t_suction * t_condensation**2 * data[8, ic] +
                            t_condensation**3 * data[9, ic] )

            ic += 1
            flow_rate[isp] = (data[0, ic] + t_suction * data[1, ic] + t_condensation * data[2, ic] +
                            t_suction**2 * data[3, ic] + t_suction * t_condensation * data[4, ic] +
                            t_condensation**2 * data[5, ic] + t_suction**3 *
                            data[6, ic] + t_suction**2 * t_condensation * data[7, ic] +
                            t_suction * t_condensation**2 * data[8, ic] +
                            t_condensation**3 * data[9, ic] )

            ic += 1
            current[isp] = (data[0, ic] + t_suction * data[1, ic] + t_condensation * data[2, ic] +
                            t_suction**2 * data[3, ic] + t_suction * t_condensation * data[4, ic] +
                            t_condensation**2 * data[5, ic] + t_suction**3 *
                            data[6, ic] + t_suction**2 * t_condensation * data[7, ic] +
                            t_suction * t_condensation**2 * data[8, ic] +
                            t_condensation**3 * data[9, ic] )


        cap = griddata(np.array(self.rpm), capacity, speed, method='linear')
        pow = griddata(np.array(self.rpm), power, speed, method='linear')
        flow = griddata(np.array(self.rpm), flow_rate, speed, method='linear')
        cur = griddata(np.array(self.rpm), current, speed, method='linear')
        coolcap = cap - pow
        cop = cap / pow


        return np.array((cap, pow, flow,cur, coolcap, cop))

    def h2 (self, Pel, h1, mdot_refC, minj, hinj):
        # Get Hotgas Enthalpy

        h1_ges = ((mdot_refC - minj) * h1 + minj * hinj) / mdot_refC

        h2_ges = h1_ges + Pel / mdot_refC

        return h2_ges




# Hier starten die Kompressoren als Objekte verf??gbar

TGH285 = TurboCor('TGH285')
TTH375 = TurboCor('TTH375')

TGH285_noEcon = TurboCor_noEcon('TGH285')
TTH375_noEcon = TurboCor_noEcon('TTH375')

AVB87DA203_50kW = PolyScroll('Belaria_100_R32_Mitsubishi_Polynome_AVB87DA203_50kW.csv', [30, 60, 100, 120])

print(AVB87DA203_50kW.calculate_fromrpm(120, -6, 35.5))
