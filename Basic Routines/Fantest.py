import numpy as np
from scipy.interpolate import griddata

FanData = np.array(((33630, 0), (29845, 80), (25915, 150), (20280, 220), (29025, 0), (25820, 61), (22440, 112),
                    (17560, 165), (23900, 0), (21260, 41), (184804, 76), (14460, 112), (18780, 0), (16075, 26),
                    (14520, 47), (11360, 69)))

FanPower = np.array((1563,
                     1938,
                     2245,
                     2550,
                     1004,
                     1255,
                     1457,
                     1663,
                     561,
                     701,
                     814,
                     929,
                     272,
                     340,
                     395,
                     451))

vol = 16000.0
dpres = 67.0

value_16T67 = griddata(FanData, FanPower, (vol, dpres), method='nearest')
value_16T67lin = griddata(FanData, FanPower, (vol, dpres), method='linear')

print(value_16T67lin, value_16T67)