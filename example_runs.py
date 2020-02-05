from airfoil import Airfoil
from foils.CSTcurve import CSTairfoil, CSTcurve


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

foil = airfoil.genAirfoil()
foil.plotAirfoil()

# 1.4 it is basically the same airfoil, but we now have plenty of options available for this generic Airfoil class object
# lets modify trailing edge to make it more realistic

foil.cutTE_XFOIL(t = .01, r = 0.25)
foil.plotAirfoil()

# several options are available for geometry modifications including: leading edge radius modification (e.g. to match required radius)
# thickenss and camber scaling (e.g. to match required maximum thickness), etc.

# 1.5 lets take a look at last part of Airfoil class- analysis
# there are options for xfoil and fluent analyses
# fluent requires software, which we might not have, so lets stick to xfoil and run our foil at re = 5e5 at angle of attack = 5 degrees:

foil.runXFOIL(alfa = 5, re=5e5)

# outputs ->

