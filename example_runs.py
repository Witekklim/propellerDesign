from airfoil import Airfoil
from foils.CSTcurve import CSTairfoil



# 1. will build the airfoil using the cst parametric definition
# 1.1. define some arbitrary parameters for cst, as vectors P and Q:
P = [0.03, .08,  .12 , .16,  .06]
Q = [.2, .25, .25,  .2, .1]
# 1.2. build airfoil of class CST, using default values of exponents: N1, N2, M1, M2 and 150 points:
airfoil = CSTairfoil(P, Q, N1 = .5, N2 = .5,  M1 = 1, M2 = 1, points = 150)
airfoil.plotAirfoil()
# 1.3. the resulting airfoil is nice and smooth, but it most likely requires some additional processing to be useful
# here i mean e.g. finite thickness trailing edge
# so lets build Airfoil class object basing on our airfoil, name it foil as airfoil is taken by CST object:
#foil = airfoil.genAirfoil()
#foil.plotAirfoil()
#foil.findCamberThickness(plot = True)

