import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import os

from Fluids import Fluid
from PumpsandFans import Pump


class Compressor:
    """ Hauptklasse """

    def __init__(self, hsuc, hdisch, hliquid, mflow):
        self.enthalpy_liquid_subcooled = hliquid
        self.enthalpy_suction = hsuc
        self.enthalpy_discharge = hdisch
        self.massflow_discharge = mflow

    def q_condenser(self):
        """ Leistung Kondensator"""
        return (self.enthalpy_discharge - self.enthalpy_liquid_subcooled) * self.massflow_discharge


class TurboCor(Compressor):
    """ Klasse TurboCor Kompressor von Danfoss
        Daten sind aus der SW von Danfoss abgeleitet. Das Original-Excel von Danfoss wurde modifiziert.
        Die Verdamfungstemperatur, die Kondensationstemperatur und die Leistung wurden mit und ohne Economizer in Schritten varriert
        und anschliessend in eine csv Datei geschrieben.
        Die Klasse TurboCor liest diese Dateien und interpoliert linear - mit und ohne Economizer und liefert das bessere Ergebnis
        self.economizer wird dann gesetzt.
    """

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
        """ Liefert """
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

        """ Bevorzugte Funktion:
            Liefert valid = get_values für alle Daten des Betriebspunkt zeropower = optimale Effizienz
        """

        valid = self.get_values(t_suction, t_condensation, self.zeropower_econ_data)

        if valid == 0:
            valid = self.get_values(t_suction, t_condensation, self.zeropower_noecon_data)
            self.economizer = 0
        else:
            self.economizer = 1
        return valid

    def calculate_frompower(self, frompower, t_suction, t_condensation):
        """ Bevorzugte Funktion:
            Liefert valid = get_values für alle Daten des Betriebspunkt frompower = prozentuale Leistung
        """

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
    """ Moduliierende Scroll Verdichter auf der Basis der Standardpolynome"""

    def __init__(self, fname: str, rpm):
        self.rpm = rpm

        # Compressor data from Polynoms with speeds
        print(os.getcwd())

        fname = "..\\Datenfiles\\" + fname
        fcsv = fname
        pddata = pd.read_csv(fcsv, sep=';')
        linedat = pddata.to_numpy(copy=True)
        self.poly_data = linedat

    def calculate_direct(self, speed, t_suction, t_condensation):
        """ Eingabe: Drehzahl, Verdampfungstemperatur, Kondensationstemperatur
             Ergebnis: nparray: Heizleistung, Kompressorleistung, Massenstrom, Stromaufnahme, Verdampferleistung, COP
             Verwendet griddata für die Interpolation"""
        data = self.poly_data
        il = 2
        co = 9

        capacity = float(data[il, co]) + data[il, co + 1] * t_suction + data[il, co + 2] * t_condensation + \
                   data[il, co + 3] * t_suction ** 2 + data[il, co + 4] * t_condensation ** 2

        cap1 = t_suction * t_condensation * (
                speed * speed * (
                data[il, co + 5] + data[il, co + 6] * t_suction + data[il, co + 7] * t_condensation) +
                speed * (
                        data[il, co + 8] + data[il, co + 9] * t_suction + data[il, co + 10] * t_condensation) +
                data[il, co + 11] + data[il, co + 12] * t_suction + data[il, co + 13] * t_condensation
        )

        cap2 = data[il, co + 14] * t_suction ** 3 + data[il, co + 15] * t_condensation ** 3

        cap3 = speed * (

                data[il, co + 16] + data[il, co + 17] * t_suction + data[il, co + 18] * t_condensation +
                data[il, co + 19] * t_suction ** 2 + data[il, co + 20] * t_condensation ** 2 +
                data[il, co + 21] * t_suction ** 3 + data[il, co + 22] * t_condensation ** 3)

        cap4 = speed ** 2 * (
                data[il, co + 23] + data[il, co + 24] * t_suction + data[il, co + 25] * t_condensation +
                data[il, co + 26] * t_suction ** 2 + data[il, co + 27] * t_condensation ** 2 +
                data[il, co + 28] * t_suction ** 3 + data[il, co + 29] * t_condensation ** 3)
        capacity = capacity + cap1 + cap2 + cap3 + cap4

        il = 1  # Achtung 1 Epower mit Drive
        co = 39


        epower = float(data[il, co]) + data[il, co + 1] * t_suction + data[il, co + 2] * t_condensation + \
                 data[il, co + 3] * t_suction ** 2 + data[il, co + 4] * t_condensation ** 2

        epower1 = t_suction * t_condensation * (
                speed * speed * (
                data[il, co + 5] + data[il, co + 6] * t_suction + data[il, co + 7] * t_condensation) +
                speed * (
                        data[il, co + 8] + data[il, co + 9] * t_suction + data[il, co + 10] * t_condensation) +
                data[il, co + 11] + data[il, co + 12] * t_suction + data[il, co + 13] * t_condensation
        )

        epower2 = data[il, co + 14] * t_suction ** 3 + data[il, co + 15] * t_condensation ** 3

        epower3 = speed * (

                data[il, co + 16] + data[il, co + 17] * t_suction + data[il, co + 18] * t_condensation +
                data[il, co + 19] * t_suction ** 2 + data[il, co + 20] * t_condensation ** 2 +
                data[il, co + 21] * t_suction ** 3 + data[il, co + 22] * t_condensation ** 3)

        epower4 = speed ** 2 * (
                data[il, co + 23] + data[il, co + 24] * t_suction + data[il, co + 25] * t_condensation +
                data[il, co + 26] * t_suction ** 2 + data[il, co + 27] * t_condensation ** 2 +
                data[il, co + 28] * t_suction ** 3 + data[il, co + 29] * t_condensation ** 3)
        epower = epower + epower1 + epower2 + epower3 + epower4

        il = 2  # Jetzt der Strom von line 2
        co = 69

        ecurrent = float(data[il, co]) + data[il, co + 1] * t_suction + data[il, co + 2] * t_condensation + \
                 data[il, co + 3] * t_suction ** 2 + data[il, co + 4] * t_condensation ** 2

        ecurrent1 = t_suction * t_condensation * (
                speed * speed * (
                data[il, co + 5] + data[il, co + 6] * t_suction + data[il, co + 7] * t_condensation) +
                speed * (
                        data[il, co + 8] + data[il, co + 9] * t_suction + data[il, co + 10] * t_condensation) +
                data[il, co + 11] + data[il, co + 12] * t_suction + data[il, co + 13] * t_condensation
        )

        ecurrent2 = data[il, co + 14] * t_suction ** 3 + data[il, co + 15] * t_condensation ** 3

        ecurrent3 = speed * (

                data[il, co + 16] + data[il, co + 17] * t_suction + data[il, co + 18] * t_condensation +
                data[il, co + 19] * t_suction ** 2 + data[il, co + 20] * t_condensation ** 2 +
                data[il, co + 21] * t_suction ** 3 + data[il, co + 22] * t_condensation ** 3)

        ecurrent4 = speed ** 2 * (
                data[il, co + 23] + data[il, co + 24] * t_suction + data[il, co + 25] * t_condensation +
                data[il, co + 26] * t_suction ** 2 + data[il, co + 27] * t_condensation ** 2 +
                data[il, co + 28] * t_suction ** 3 + data[il, co + 29] * t_condensation ** 3)
        ecurrent = ecurrent + ecurrent1 + ecurrent2 + ecurrent3 + ecurrent4

        il = 2  # Jetzt die Kompressorleistung von line 2 für die Heizleistung
        co = 39

        comppower = float(data[il, co]) + data[il, co + 1] * t_suction + data[il, co + 2] * t_condensation + \
                 data[il, co + 3] * t_suction ** 2 + data[il, co + 4] * t_condensation ** 2

        comppower1 = t_suction * t_condensation * (
                speed * speed * (
                data[il, co + 5] + data[il, co + 6] * t_suction + data[il, co + 7] * t_condensation) +
                speed * (
                        data[il, co + 8] + data[il, co + 9] * t_suction + data[il, co + 10] * t_condensation) +
                data[il, co + 11] + data[il, co + 12] * t_suction + data[il, co + 13] * t_condensation
        )

        comppower2 = data[il, co + 14] * t_suction ** 3 + data[il, co + 15] * t_condensation ** 3

        comppower3 = speed * (

                data[il, co + 16] + data[il, co + 17] * t_suction + data[il, co + 18] * t_condensation +
                data[il, co + 19] * t_suction ** 2 + data[il, co + 20] * t_condensation ** 2 +
                data[il, co + 21] * t_suction ** 3 + data[il, co + 22] * t_condensation ** 3)

        comppower4 = speed ** 2 * (
                data[il, co + 23] + data[il, co + 24] * t_suction + data[il, co + 25] * t_condensation +
                data[il, co + 26] * t_suction ** 2 + data[il, co + 27] * t_condensation ** 2 +
                data[il, co + 28] * t_suction ** 3 + data[il, co + 29] * t_condensation ** 3)
        comppower = comppower + comppower1 + comppower2 + comppower3 + comppower4

        il = 2  # Jetzt der Massenstrom von line 2
        co = 99

        massflow = float(data[il, co]) + data[il, co + 1] * t_suction + data[il, co + 2] * t_condensation + \
                 data[il, co + 3] * t_suction ** 2 + data[il, co + 4] * t_condensation ** 2

        massflow1 = t_suction * t_condensation * (
                speed * speed * (
                data[il, co + 5] + data[il, co + 6] * t_suction + data[il, co + 7] * t_condensation) +
                speed * (
                        data[il, co + 8] + data[il, co + 9] * t_suction + data[il, co + 10] * t_condensation) +
                data[il, co + 11] + data[il, co + 12] * t_suction + data[il, co + 13] * t_condensation
        )

        massflow2 = data[il, co + 14] * t_suction ** 3 + data[il, co + 15] * t_condensation ** 3

        massflow3 = speed * (

                data[il, co + 16] + data[il, co + 17] * t_suction + data[il, co + 18] * t_condensation +
                data[il, co + 19] * t_suction ** 2 + data[il, co + 20] * t_condensation ** 2 +
                data[il, co + 21] * t_suction ** 3 + data[il, co + 22] * t_condensation ** 3)

        massflow4 = speed ** 2 * (
                data[il, co + 23] + data[il, co + 24] * t_suction + data[il, co + 25] * t_condensation +
                data[il, co + 26] * t_suction ** 2 + data[il, co + 27] * t_condensation ** 2 +
                data[il, co + 28] * t_suction ** 3 + data[il, co + 29] * t_condensation ** 3)
        massflow = massflow + massflow1 + massflow2 + massflow3 + massflow4

        il = 2  # Jetzt die Discharge Temperature von line 2
        co = 129

        Tdischarge = float(data[il, co]) + data[il, co + 1] * t_suction + data[il, co + 2] * t_condensation + \
                 data[il, co + 3] * t_suction ** 2 + data[il, co + 4] * t_condensation ** 2

        Tdischarge1 = t_suction * t_condensation * (
                speed * speed * (
                data[il, co + 5] + data[il, co + 6] * t_suction + data[il, co + 7] * t_condensation) +
                speed * (
                        data[il, co + 8] + data[il, co + 9] * t_suction + data[il, co + 10] * t_condensation) +
                data[il, co + 11] + data[il, co + 12] * t_suction + data[il, co + 13] * t_condensation
        )

        Tdischarge2 = data[il, co + 14] * t_suction ** 3 + data[il, co + 15] * t_condensation ** 3

        Tdischarge3 = speed * (

                data[il, co + 16] + data[il, co + 17] * t_suction + data[il, co + 18] * t_condensation +
                data[il, co + 19] * t_suction ** 2 + data[il, co + 20] * t_condensation ** 2 +
                data[il, co + 21] * t_suction ** 3 + data[il, co + 22] * t_condensation ** 3)

        Tdischarge4 = speed ** 2 * (
                float(data[il, co + 23]) + data[il, co + 24] * t_suction + data[il, co + 25] * t_condensation +
                data[il, co + 26] * t_suction ** 2 + data[il, co + 27] * t_condensation ** 2 +
                data[il, co + 28] * t_suction ** 3 + data[il, co + 29] * t_condensation ** 3)
        Tdischarge = Tdischarge + Tdischarge1 + Tdischarge2 + Tdischarge3 + Tdischarge4

        qheat = capacity + comppower
        cop = qheat / epower

        return np.array((qheat, epower, massflow, ecurrent, capacity, cop, Tdischarge))

    def calculate_fromrpm(self, speed, t_suction, t_condensation):
        """ Eingabe: Drehzahl, Verdampfungstemperatur, Kondensationstemperatur
            Ergebnis: nparray: Heizleistung, Kompressorleistung, Massenstrom, Stromaufnahme, Verdampferleistung, COP
            Verwendet griddata für die Interpolation"""
        data = self.poly_data
        ishift = 7
        rpm = self.rpm
        capacity = np.empty(4, float)  # nur fuer die Dimension der Liste
        power = np.empty(4, float)  # nur fuer die Dimension der Liste
        flow_rate = np.empty(4, float)  # nur fuer die Dimension der Liste
        current = np.empty(4, float)  # nur fuer die Dimension der Liste

        for isp in range(0, len(self.rpm)):
            ic = ishift * isp + 1
            capacity[isp] = (data[0, ic] + t_suction * data[1, ic] + t_condensation * data[2, ic] +
                             t_suction ** 2 * data[3, ic] + t_suction * t_condensation * data[4, ic] +
                             t_condensation ** 2 * data[5, ic] + t_suction ** 3 *
                             data[6, ic] + t_suction ** 2 * t_condensation * data[7, ic] +
                             t_suction * t_condensation ** 2 * data[8, ic] +
                             t_condensation ** 3 * data[9, ic])
            ic += 1
            power[isp] = (data[0, ic] + t_suction * data[1, ic] + t_condensation * data[2, ic] +
                          t_suction ** 2 * data[3, ic] + t_suction * t_condensation * data[4, ic] +
                          t_condensation ** 2 * data[5, ic] + t_suction ** 3 *
                          data[6, ic] + t_suction ** 2 * t_condensation * data[7, ic] +
                          t_suction * t_condensation ** 2 * data[8, ic] +
                          t_condensation ** 3 * data[9, ic])

            ic += 1
            flow_rate[isp] = (data[0, ic] + t_suction * data[1, ic] + t_condensation * data[2, ic] +
                              t_suction ** 2 * data[3, ic] + t_suction * t_condensation * data[4, ic] +
                              t_condensation ** 2 * data[5, ic] + t_suction ** 3 *
                              data[6, ic] + t_suction ** 2 * t_condensation * data[7, ic] +
                              t_suction * t_condensation ** 2 * data[8, ic] +
                              t_condensation ** 3 * data[9, ic])

            ic += 1
            current[isp] = (data[0, ic] + t_suction * data[1, ic] + t_condensation * data[2, ic] +
                            t_suction ** 2 * data[3, ic] + t_suction * t_condensation * data[4, ic] +
                            t_condensation ** 2 * data[5, ic] + t_suction ** 3 *
                            data[6, ic] + t_suction ** 2 * t_condensation * data[7, ic] +
                            t_suction * t_condensation ** 2 * data[8, ic] +
                            t_condensation ** 3 * data[9, ic])

        cap = griddata(np.array(self.rpm), capacity, speed, method='linear')
        pow = griddata(np.array(self.rpm), power, speed, method='linear')
        flow = griddata(np.array(self.rpm), flow_rate, speed, method='linear')
        cur = griddata(np.array(self.rpm), current, speed, method='linear')
        coolcap = cap - pow
        cop = cap / pow

        return np.array((cap, pow, flow, cur, coolcap, cop))

    def h2(self, Pel, h1, mdot_refC, minj, hinj):
        # Get Hotgas Enthalpy

        h1_ges = ((mdot_refC - minj) * h1 + minj * hinj) / mdot_refC

        h2_ges = h1_ges + Pel / mdot_refC

        return h2_ges

    def h2a(self, Pel, h1, mdot_refC, minj, hinj):
        # Get Hotgas Enthalpy

        h1_ges = (mdot_refC * h1 + minj * hinj) / (mdot_refC + minj)

        h2_ges = h1_ges + Pel / (mdot_refC + minj)

        return h2_ges


class fix_PolyScroll(PolyScroll):
    """ Fixed speed Scroll Verdichter mit Polynom"""

    def __init__(self, fname: str, rpm):
        self.rpm = rpm

        # Compressor data from Polynoms with speeds
        fname = "..\\Datenfiles\\" + fname
        fcsv = fname
        pddata = pd.read_csv(fcsv, sep=';')
        linedat = pddata.to_numpy(copy=True)
        self.poly_data = linedat

    def calculate_fromrpm(self, speed, t_suction, t_condensation):
        """ Eingabe: Drehzahl (dummy), Verdampfungstemperatur, Kondensationstemperatur
             Ergebnis: nparray: Heizleistung, Kompressorleistung, Massenstrom, Stromaufnahme, Verdampferleistung, COP, Heissgastemperatur (sofern vorhanden)
            """

        data = self.poly_data
        ishift = 0
        rpm = self.rpm
        capacity = np.empty(4, float)  # nur fuer die Dimension der Liste
        power = np.empty(4, float)  # nur fuer die Dimension der Liste
        flow_rate = np.empty(4, float)  # nur fuer die Dimension der Liste
        current = np.empty(4, float)  # nur fuer die Dimension der Liste

        cap = (data[0, 1] + t_suction * data[1, 1] + t_condensation * data[2, 1] +
               t_suction ** 2 * data[3, 1] + t_suction * t_condensation * data[4, 1] +
               t_condensation ** 2 * data[5, 1] + t_suction ** 3 *
               data[6, 1] + t_suction ** 2 * t_condensation * data[7, 1] +
               t_suction * t_condensation ** 2 * data[8, 1] +
               t_condensation ** 3 * data[9, 1])

        power = (data[0, 2] + t_suction * data[1, 2] + t_condensation * data[2, 2] +
                 t_suction ** 2 * data[3, 2] + t_suction * t_condensation * data[4, 2] +
                 t_condensation ** 2 * data[5, 2] + t_suction ** 3 *
                 data[6, 2] + t_suction ** 2 * t_condensation * data[7, 2] +
                 t_suction * t_condensation ** 2 * data[8, 2] +
                 t_condensation ** 3 * data[9, 2])

        flow = (data[0, 3] + t_suction * data[1, 3] + t_condensation * data[2, 3] +
                t_suction ** 2 * data[3, 3] + t_suction * t_condensation * data[4, 3] +
                t_condensation ** 2 * data[5, 3] + t_suction ** 3 *
                data[6, 3] + t_suction ** 2 * t_condensation * data[7, 3] +
                t_suction * t_condensation ** 2 * data[8, 3] +
                t_condensation ** 3 * data[9, 3])

        cur = (data[0, 4] + t_suction * data[1, 4] + t_condensation * data[2, 4] +
               t_suction ** 2 * data[3, 4] + t_suction * t_condensation * data[4, 4] +
               t_condensation ** 2 * data[5, 4] + t_suction ** 3 *
               data[6, 4] + t_suction ** 2 * t_condensation * data[7, 4] +
               t_suction * t_condensation ** 2 * data[8, 4] +
               t_condensation ** 3 * data[9, 4])

        tdis = (data[0, 5] + t_suction * data[1, 5] + t_condensation * data[2, 5] +
                t_suction ** 2 * data[3, 5] + t_suction * t_condensation * data[4, 5] +
                t_condensation ** 2 * data[5, 5] + t_suction ** 3 *
                data[6, 5] + t_suction ** 2 * t_condensation * data[7, 5] +
                t_suction * t_condensation ** 2 * data[8, 5] +
                t_condensation ** 3 * data[9, 5])

        coolcap = cap - power
        cop = cap / power

        return np.array((cap, power, flow, cur, coolcap, cop, tdis))


# Hier starten die Kompressoren als Objekte verfügbar


VZN175 = PolyScroll('Belaria_65_R290_Danfoss_Polynome_VZN175_69kW.csv', [30, 140])
# VZN175 = PolyScroll('qwertz.csv', [30, 140])
# s = VZN175.calculate_direct(120, -16.0, 56.0)
# print(s)
# s = VZN175.calculate_direct(120, -16.0, 57.0)
# print(s)
# s = VZN175.calculate_direct(120, -1640, 57.0)
# print(s)
# s = VZN175.calculate_direct(120, -7.0, 37.0)
# print(s)
# s = VZN175.calculate_direct(120, 2.0, 37.0)
# print(s)


TGH285 = TurboCor('TGH285')
"""TGH-285 mit R134a/R515A"""
TTH375 = TurboCor('TTH375')
"""TTH-375 mit T1234ze/513B"""

TGH285_noEcon = TurboCor_noEcon('TGH285')
"""TGH-285 mit R134a/R515A ohne Economizer"""
TTH375_noEcon = TurboCor_noEcon('TTH375')
"""TTH-375 mit T1234ze/513B ohne Economizer"""

AVB87DA203_50kW = PolyScroll('Belaria_100_R32_Mitsubishi_Polynome_AVB87DA203_50kW.csv', [30, 60, 100, 120])
""" Mitsubishi Kompressor für Belaria (100) mit R32 - nicht mehr verwendet - modulierend Drehzahlen [30, 60, 100, 120]"""
DSG480_4_120kW = fix_PolyScroll('R515B_System_Polynome_DSG480_4_120kW.csv', 58)
""" Fixed Speed Danfoss DSG480 - preliminary data R515B"""
DSG480_4R_120kW = fix_PolyScroll('R1234ze_System_Polynome_DSG480_4_120kW.csv', 58)
""" Fixed Speed Danfoss DSG480 - preliminary data R1234ze"""
YH33K1G_33kW = fix_PolyScroll('R290_Emerson_YH33K1G.csv', 50)
""" Emerson fixed speed für Dual Twin Konzept Belaria 100 in R290"""

print(AVB87DA203_50kW.calculate_fromrpm(120, -6, 35.5))
print(DSG480_4R_120kW.calculate_fromrpm(58, -6, 35.5))
print(DSG480_4_120kW.calculate_fromrpm(58, -6, 35.5))
print(YH33K1G_33kW.calculate_fromrpm(50, -7, 35))
print(VZN175.calculate_direct(120, 2.0, 37.0))