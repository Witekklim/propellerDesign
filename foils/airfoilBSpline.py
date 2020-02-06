# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 14:30:14 2019

@author: WK5521

knotsU = np.array([[0,0], [0., .04], [0.2, .1], [.4, .12],[.6,.1], [.8, .05] ,[1,0] ])
knotsL = np.array([[0,0], [0., -.025], [0.2, -.05], [.4, -0.06],[.6, -.01],[.8, .01], [1,0] ])
foil = BSplineAirfoil(knotsU, knotsL)
foil.plotAirfoil(True)
foil.saveTex('LB')


a = foil.genAirfoil()
a.plotAirfoil()
a.findCamberThickness(True)

a.runFluent(2,.2,1,'bspline1')


-------------------------------------------------------------------------
BS airofil from 10 parameters: 
x = np.asarray([ 0.02774923,  0.09707878,  0.08142534,  0.06446338,  0.02688518,-0.0535    , -0.12202614, -0.11585406, -0.08      , -0.03      ])
knotsU = np.asarray([[0, 0], [0, x[0]], [0.2, x[1]], [0.4, x[2]], [.6, x[3]], [.8, x[4]], [1, 0]])
knotsL = np.asarray([[0, 0], [0, x[0] +x[5]], [0.2,x[1] +  x[6]], [0.4,x[2] +  x[7]], [.6,x[3] + x[8]], [.8,x[4] +  x[9]], [1, 0]])
foil = BSplineAirfoil(knotsU, knotsL)
foil.plotAirfoil(True)
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import sys
import os

sys.path.append(os.getcwd().strip('\\foils'))
#sys.path.append("E:\propeller\python")

from airfoil import Airfoil



class BSplineAirfoil(object):
    """  build a bspline based airfoil basing on two sets of kntos: for upper and lower surfaces """
    def __init__(self, knotsU, knotsL):
        self.knotsU = knotsU
        self.knotsL = knotsL
        
        self.curveU = curve(self.knotsU, 100)
        self.curveL = curve(self.knotsL, 100)
        
        self.xU = self.curveU.curve[0]
        self.yU = self.curveU.curve[1]
        self.xL = self.curveL.curve[0]
        self.yL = self.curveL.curve[1]
        
        self.x = np.concatenate ((np.flip(self.xL, axis = 0)[:-1], self.xU ))
        self.y = np.concatenate ((np.flip(self.yL, axis = 0)[:-1], self.yU )) 
        
#        plt.plot(self.x, self.y)
        
        self.fileFig = r'C:\Users\wk5521\Documents\python\saved_plots\airfoil'
        
    def plotAirfoil(self, saveFig = False):
        
        plt.figure(figsize = (6, 3), dpi = 200)
        plt.plot(self.xU, self.yU, 'k-', label = 'airfoil')
        plt.plot(self.xL, self.yL, 'k-')
        plt.plot(self.knotsU[:,0], self.knotsU[:,1], 'ko--', mfc='none', linewidth = 0.7, label = 'control polygon')
        plt.plot(self.knotsL[:,0], self.knotsL[:,1], 'ko--', mfc='none', linewidth = 0.7)
        plt.legend()
        plt.xlabel('$x/c$')
        plt.ylabel('$y/c$')
        plt.tight_layout()
        if saveFig:
            plt.savefig(self.fileFig, dpi = 1000)
        plt.show()
        
    def saveTex(self,name):
        X = np.append(self.knotsU.reshape(-1,2),np.flip( self.knotsL[:-1].reshape(-1,2), axis = 0), axis = 0)
        np.savetxt(r'E:\propeller\python\wing3d\tex-plots\knots{}.txt'.format(name), X)
        
        
        X = np.append(self.xU.reshape(-1,1), self.yU.reshape(-1,1), axis = 1)
        y = np.append(self.xL.reshape(-1,1), self.yL.reshape(-1,1), axis = 1)
        X = np.append(X, np.flip(y, axis = 0), axis = 0)
        np.savetxt(r'E:\propeller\python\wing3d\tex-plots\points{}.txt'.format(name), X)
        
        
        
        
    def genAirfoil(self, t = 0.005):
#        self.(self, ftype = 'ICEM', filein = None, x = None, y = None,  T_req = None, camber = None,
#                 chord = None, beta = None, z = 0, fileoutICEM = None, t = 0, dx = 0, dy = 0, split = False, origin = 0,
#                 verbose = False)
        
        return Airfoil('XY', x = self.x , y=self.y, t = t)
        
        
        
        
        
        
        
class curve:
    """
    cur = curve( np.array([[0,1],[1,2],[1,4]]))

    self.curve defines bspline from input knots
    
    """
    def __init__(self, knots, points = 100):
        """ finds bspline curve for specified knots
        inputs: knots coords as numpy array  [  [x0,y0],[x1,y1],...,[xn-1,yn-1]  ]
        attributes:
            (x,y):      knots used for spline definition 
            curve:      (x,y) coordinates of resulting curve (size specified by points)
                        [ [x0,x1,x2,...,xn-1], [y0,y1,y2,...,yn-1] ]
        """
        self.x=knots[:,0]
        self.y=knots[:,1]
                
        l=len(self.x)  
        t=np.linspace(0,1,l-2,endpoint=True)
        t=np.append([0,0,0],t)
        t=np.append(t,[1,1,1])
        tck=[t,[self.x,self.y], 3]
        u3=np.linspace(0,1,(max(l*2,points)),endpoint=True)
        self.curve = interpolate.splev(u3, tck, 0)
        
        
     
        
  
