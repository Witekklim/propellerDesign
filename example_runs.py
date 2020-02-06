"""
following file contains example commands used to use the library
it is focused on airfoil design, which uses parameteric definitions of airfoils and runs optimization procedures
"""
# 1. will build the airfoil using the cst parametric definition
# 1.1. define some arbitrary parameters for cst, as vectors P and Q:
P = [.03, .08, .12, .16, .06]
Q = [ .2, .25, .25,  .2, .1]

# here P is camber, Q is thickness distribution. more parameters can be added to allow more geometry flexibility
# 1.2. build airfoil of class CST, using default values of exponents: N1, N2, M1, M2 and 150 points:
from foils.airfoilCST import CSTairfoil, CSTcurve
airfoil = CSTairfoil(P, Q, N1 = .5, N2 = .5,  M1 = 1, M2 = 1, points = 150)
airfoil.plotAirfoil()

# 1.3. the resulting airfoil is nice and smooth, but it most likely requires some additional processing to be useful
# here i mean e.g. finite thickness trailing edge
# so lets build Airfoil class object basing on our airfoil, name it foil as airfoil is taken by CST object:
foil = airfoil.genAirfoil()
foil.plotAirfoil()

# 1.4. it is basically the same airfoil, but we now have plenty of options available for this generic Airfoil class object
# lets modify trailing edge to make it more realistic
foil.cutTE_XFOIL(t = .01, r = 0.25)
foil.plotAirfoil()

# several options are available for geometry modifications including: leading edge radius modification (e.g. to match required radius)
# thickenss and camber scaling (e.g. to match required maximum thickness), etc.

# 1.5. lets take a look at last part of Airfoil class- analysis
# there are options for xfoil and fluent analyses
# fluent requires software, which we might not have, so lets stick to xfoil and run our foil at re = 5e5 at angle of attack = 5 degrees:
foil.runXFOIL(alfa = 5, re=5e5)

# outputs -> (5.0 1.0132 0.01035 -0.1013), which means 1.0132 lift coefficient, 0.01035 drag coefficient and -0.1013 pitching moment coefficient
# a several options for running xfoil are defined as well, which allow to quickly run airfoil polar or run for desired lift coefficient

#==================================================================
# 2. similar procedure can be used to build b-spline based airfoil
from foils.airfoilBSpline import BSplineAirfoil
import numpy as np
# 2.1. lets define vectors of knots which define splines
knotsU = np.array([[0,0], [0., .04], [0.2, .1], [.4, .12],[.6,.1], [.8, .05] ,[1,0] ])
knotsL = np.array([[0,0], [0., -.025], [0.2, -.05], [.4, -0.06],[.6, -.01],[.8, .01], [1,0] ])

# and build BSplineAirfoil class foil:
airfoilb = BSplineAirfoil(knotsU, knotsL)
airfoilb.plotAirfoil()

# plot presents control polygon and resulting airfoil
# 2.2. the same procedure as for cst foil can be used to bspline foil:
foilb = airfoilb.genAirfoil()
foilb.plotAirfoil(name = 'bspline foil as Airfoil object')
foilb.cutTE_XFOIL(t = .01, r = 0.25)
foilb.plotAirfoil(name = 'thickened te of bspline foil')

# needless to mention, this library allows processing of airfoil to match required geoemtry parameters of any form: 
# chord, maximum thickness, leading edge radius, trailing edge thickness, camber and thickness scaling, rotating, etc.
# hence it can be used to define airfoil stacking in 3d wing-type design cases.

#==================================================================
# 3. run airfoil design using evolutionary process
from foils.airfoilDesigner import airfoilDesigner
xs = [.02, .3, .9]
weights = [1, 1, 1]
maxiter = 80

design = airfoilDesigner([1],'CST', 10, ts = [.04, .11, .015], xs = xs, weights = weights)
design.runOptimization(maxiter = maxiter)
