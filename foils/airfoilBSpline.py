# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 14:30:14 2019
@author: Witold Klimczyk

file contains BSplineAirfoil class, which builds bspline airofils from knots for upper nad lower curves
based on class curve, which builds bsplines with scipy implemnetation

"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import sys
import os

sys.path.append(os.getcwd().strip('\\foils'))
#sys.path.append("E:\propeller\python")

from airfoil import Airfoil



class BSplineAirfoil():
    """  build a bspline based airfoil basing on two sets of kntos: for upper and lower surfaces """
    def __init__(self, knotsU, knotsL, fileFig = None):
        self.knotsU = knotsU
        self.knotsL = knotsL
        self.fileFig = fileFig

        self.curveU = curve(self.knotsU, 100)
        self.curveL = curve(self.knotsL, 100)
        
        self.xU = self.curveU.curve[0]
        self.yU = self.curveU.curve[1]
        self.xL = self.curveL.curve[0]
        self.yL = self.curveL.curve[1]
        
        self.x = np.concatenate ((np.flip(self.xL, axis = 0)[:-1], self.xU ))
        self.y = np.concatenate ((np.flip(self.yL, axis = 0)[:-1], self.yU )) 
        
    def plotAirfoil(self, name = 'b-spline airfoil', saveFig = False):
        
        plt.figure(figsize = (6, 3), dpi = 200)
        plt.plot(self.xU, self.yU, 'k-', label = 'airfoil')
        plt.plot(self.xL, self.yL, 'k-')
        plt.plot(self.knotsU[:,0], self.knotsU[:,1], 'ko--', mfc='none', linewidth = 0.7, label = 'control polygon')
        plt.plot(self.knotsL[:,0], self.knotsL[:,1], 'ko--', mfc='none', linewidth = 0.7)
        plt.legend()
        plt.xlabel('$x/c$')
        plt.ylabel('$y/c$')
        plt.title(name)
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
        
    def genAirfoil(self, t = 0):
        return Airfoil('XY', x = self.x , y=self.y, t = t)
        
        
        
        
        
        
        
class curve():
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
        
        
     
        
  
