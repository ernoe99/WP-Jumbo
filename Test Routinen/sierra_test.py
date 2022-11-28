import math

import numpy as np
import pandas as pd
import os
from scipy.interpolate import griddata
import Fluids

import HeatExchanger
fname = "CoilStudioAuslegungVerdampfer4rowfuerFP063_SH5"

dname = "..\\Datenfiles\\" + fname  # Vielleicht hier ohne "..\\"
fcsv = dname + '.csv'
fVcsv = dname + 'V.csv'  # Verdampferdaten
print(os.getcwd())
pddata = pd.read_csv(fcsv, sep=';')
linedat = pddata.to_numpy(copy=True)
CondData = linedat[3:, 1:5]    # Schnittstelle 4 Kondensationstemp. Verdampfungstemp. Lufttemp.
tcdata = linedat[3:, 1:2]
tevapdata = linedat[3:, 2:3]
tairdata = linedat[3:, 3:4]

airvolData = linedat[3:, 23:24]   # Volume airflow
airvelData = linedat[3:, 24:25]   # Volume airflow
superheatData = linedat[3:, 37:38]   # Volume airflow
subcoolData = linedat[3:, 39:40]   # Volume airflow
CondPowerSensible = linedat[3:, 47:48]   # Leistung trocken - Spalte 49
Condpressuredrop_Air = linedat[3:, 24:25]   # Druckverlust Luft - Spalte 27
Condpressuredrop_Ref = linedat[3:, 40:41]   # Druckverlust Refrigerant - Spalte 43
Condcharge_Ref = linedat[3:, 41:42]   # Inhalt Refrigerant - Spalte 44
Condflow_Ref = linedat[3:, 39:40]   # Massenstrom Refrigerant - Spalte 42

evapData = linedat[3:, 1:5]    # Schnittstelle 4 Kondensationstemp. Verdampfungstemp. Lufttemp., Luftvolumen, Superheat, Subcooling

Tevap = -16
Tcond = 56
Tair = -7
Airvolume = 24600.0
Superheat = 0
Subcool = 1


power = float(griddata(evapData, CondPowerSensible, (Tcond, Tevap,  Tair, Airvolume), method='linear'))

print(power)

Subcool = 20

power = float(griddata(evapData, CondPowerSensible, (Tcond, Tevap, Tair, Airvolume), method='linear'))

print(power)


fname = "CoilStudioAuslegungVerdampfer4rowfuerFP063_SH0SC"

dname = "..\\Datenfiles\\" + fname
fcsv = dname + '.csv'
fVcsv = dname + 'V.csv'  # Verdampferdaten
print(os.getcwd())
pddata = pd.read_csv(fcsv, sep=';')
linedat = pddata.to_numpy(copy=True)
CondData = linedat[3:, 1:6]    # Schnittstelle 4 Kondensationstemp. Verdampfungstemp. Lufttemp.
tcdata = linedat[3:, 1:2]
tevapdata = linedat[3:, 2:3]
tairdata = linedat[3:, 3:4]

airvolData = linedat[3:, 23:24]   # Volume airflow
airvelData = linedat[3:, 24:25]   # Volume airflow
superheatData = linedat[3:, 37:38]   # Volume airflow
subcoolData = linedat[3:, 39:40]   # Volume airflow
CondPowerSensible = linedat[3:, 47:48]   # Leistung trocken - Spalte 49
Condpressuredrop_Air = linedat[3:, 24:25]   # Druckverlust Luft - Spalte 27
Condpressuredrop_Ref = linedat[3:, 40:41]   # Druckverlust Refrigerant - Spalte 43
Condcharge_Ref = linedat[3:, 41:42]   # Inhalt Refrigerant - Spalte 44
Condflow_Ref = linedat[3:, 39:40]   # Massenstrom Refrigerant - Spalte 42

evapData = linedat[3:, 1:6]    # Schnittstelle 4 Kondensationstemp. Verdampfungstemp. Lufttemp., Luftvolumen, Superheat, Subcooling

Tevap = -16
Tcond = 56
Tair = -7
Airvolume = 24600.0
Superheat = 0
Subcool = 1


power = float(griddata(evapData, CondPowerSensible, (Tcond, Tevap,  Tair, Airvolume, Subcool), method='linear'))

print(power)

Subcool = 20

power = float(griddata(evapData, CondPowerSensible, (Tcond, Tevap, Tair, Airvolume, Subcool), method='linear'))

print(power)


