from airfoil import airfoil
from foils import CSTcurve.CSTairfoil

P = [0.03, .08,  .12 , .16,  .06]
Q = [.16, .18, .2,  .11, .1]

airfoil = CSTairfoil(P, Q, N1 = .5, N2 = .5,  M1 = 1, M2 = 1, points = 150)
airfoil.plotAirfoil()
