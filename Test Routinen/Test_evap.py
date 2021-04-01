from Fluids import Brine
from HeatExchanger import Evaporator_fitted


CSF2880_fitted = Evaporator_fitted(306, 0.4183, 250, (-13.3, 24.01, -0.2405))
b1 = Brine(0, 20)

pevap = -200.0

CSF2880_fitted.solvebalance(b1, pevap)
