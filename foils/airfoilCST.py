# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 08:45:46 2019

@author: Witold Klimczyk
N for leading edge, M for trailing edge
N1: derivative -> 0;  -> 0.5 -> infty, >0.5 shit, keep it below 1, but > 0.7
N2: want smooth LE -> N2 = 0.5, want sharp LE -> N2 = 1
M1: want derivative 0 -> M1 = 1, icreasae derivative by setting M1 -> 0, direct derivative control with M1 = 0
M2: keep at 0.5.

"""
import numpy as np
import matplotlib.pyplot as plt
import os

import sys
path = os.getcwd()
sys.path.append(os.getcwd().strip('\\foils'))
#sys.path.append("E:\propeller\python")
from airfoil import Airfoil

from subprocess import run, PIPE

        
class CSTairfoil():
    def __init__(self, P, Q, N1 = 0.5 , N2 = 0.5, M1 = 1, M2 = 1, points = 160, N = None, verbose = False, option = 2, workingdir = None, xfoildir = None):
        """
        INPUTS:
        P and Q:     are vectors of variables for top and bottom airfoil curves
        N and M:     are parameters for defining smoothness of LE and TE
        points:      is number of points for each curve
        N:           is used for Panukl 3d wing generation, where N stands for section ID.
        option:      specifies whether foil top and bottom curves are defined explicitly or via camber + thickness 
        
        attibutes:
        top:         cst class curve defining top curve     attributes: top.x, top.y define coords
        bottom:      cst class curve defining bottom curve  
        """
        self.maxthickness = None
        self.area = None
        
        self.verbose = verbose
        
        self.name = 'P,'
        for i in range(len(P)):
            P[i] = round(P[i], 3)
            self.name += '{},'.format(P[i])
        self.name += 'Q,'
        for i in range(len(Q)-1):
            Q[i] = round(Q[i], 3)
            self.name += '{},'.format(Q[i])
        self.name += '{}'.format(P[-1])
            
        
        self.xfoilpath = 'E:/propeller/XFOIL/xfoil.exe'
        self.workingdir = workingdir if workingdir != None else os.getcwd()
        self.xfoildir = xfoildir if xfoildir != None else os.getcwd().strip('\\python')
        self.xfoilthickness = self.xfoildir + r'\thicknessfoil.txt'     # foil with thickness distribution
        self.xfoilcamber = self.xfoildir + r'\camberfoil.txt'         # foil with camber distribution
        self.cambertxt = 'camberdist.txt'                               # file with camber distribution
        self.xfoil = self.xfoildir + r'\cstfoil.txt'                  # final foil
        

        # build foil
        if option == 1:
            self.top = CSTcurve(P,N1,M1, points)
            self.bottom = CSTcurve(Q,N2,M2, points)
        elif option == 2:
            self.camber = CSTcurve(P, N1, M1, points)
            self.thick = CSTcurve(Q, N2, M2, points)
            
            topy = self.camber.y + np.array(self.thick.y)#/2
            boty = self.camber.y - np.array(self.thick.y)#/2
            
            self.top = CSTcurve(P, N1, M1, points)
            self.bottom = CSTcurve(Q, N2, M2, points)
            self.top.y = topy
            self.bottom.y = boty
        else:
            print('specify option for airfoil build:\noption = 1 for top and bottom separately\noption = 2 for camber plus thickness definition')
        
        self.saveXFOIL()
        
    
    def plotAirfoil(self):
        
        fig = plt.figure(dpi = 200)

        plt.plot(self.top.x, self.top.y, 'k-', label = 'top', linewidth = 1)
        plt.plot(self.bottom.x, self.bottom.y, 'k-', label = 'top', linewidth = 1)
        plt.plot(self.camber.x, self.camber.y, 'k--', label = 'camber', linewidth = 1)
        plt.xlabel('$x/c$')
        plt.ylabel('$y/c$')
        plt.axis('equal')
        plt.grid('major', linewidth = .2)
        plt.tight_layout()
        plt.show()
        
    def genAirfoil(self, t=0):
        """ generate airfoil"""
        return Airfoil('XFOIL', self.xfoil, t=t)
        
    def saveXFOIL(self, file=None):
        """  saves airfoil coords to .txt file with specified path   """
        # close trailing edge
        xfoilbottom = self.bottom.y
        xfoilbottom[-1] = 0
        xfoiltop = self.top.y
        xfoiltop[-1] = 0
        
        # save coords to file
        if file is None:
            file= self.xfoil
        with open(file, "w") as text_file:
            print("airfoilCST", file = text_file)
            for i in range(len(self.bottom.x)-1,0,-1):
                print("  {} {}".format(self.bottom.x[i], xfoilbottom[i]), file=text_file)
            for i in range(len(self.top.x)):
                print("  {} {}".format(self.top.x[i], xfoiltop[i]), file=text_file)
                
    def saveXFOILThickness(self, file=None):
        """  saves airfoil coords to .txt file with specified path   """
        # close trailing edge
        xfoilbottom = self.bottom.y
        xfoilbottom[-1] = 0
        xfoiltop = self.top.y
        xfoiltop[-1] = 0
        N_top = len(self.top.x)
        N_bot = len(self.bottom.x)
        
        # save coords to file
        if file is None:
            file = self.xfoilthickness
        with open(file, "w") as text_file:
            print("airfoilCST", file = text_file)
            for i in range(len(self.bottom.x)-1,0,-1):
                print("  {} {}".format(self.thick.x[i], -self.thick.y[i]), file=text_file)
            for i in range(len(self.top.x)):
                print("  {} {}".format(self.thick.x[i], self.thick.y[i]), file=text_file)
                
                
    def addCamber_XFOIL(self):
        """ modifies airfoil using xfoil to maintain camber

            t: thickness
            r: blending radius
        """
        # save foil with thickness distribution
        self.saveXFOILThickness()
        
        # save xfoil file with camber definition
        self.saveXFOIL(self.xfoilcamber)
        
        command = 'load ' + self.xfoilcamber + '\npane\ngdes\ncamb\nwrtc\n{}\ny\n'.format(self.cambertxt) + '\n\n\nquit\n'
        p = run([self.xfoilpath], stdout=PIPE,  input=command, encoding='ascii', shell = False)

        # load thickness foil and apply camber to it, save under final foil name
        command1 = 'load {}\npane\ngdes\ncamb\nrdac\n{}\nadd\n\n\npane\nsave {}\ny\nquit\n'.format( self.xfoilthickness, self.cambertxt, self.xfoil)
        p = run([self.xfoilpath], stdout=PIPE,  input=command1, encoding='ascii', shell = False)
        
        if self.verbose:
            print('succesfully constructed foil')        
                




class CSTcurve():
    """
    evaluates bernstein polynomials to form a curve for given parameters

    INPUTS:
    P:          vector of parameters of any length e.g. P = [1,1,1,1,1]
    N1, N2:     le and te exponents for derivative control, default set (0.5,1)
    
    attributes:
    x,y:        coordinates of curve
    """
    
    def __init__(self,P,N1,N2, points):
        self.P = P
        self.N1 = N1
        self.N2 = N2
        self.x = []
        self.y = []
        self.CST_curve(P, points)
        
    def Bernstein(self,i,p,phi):
    # returns value of bernstein for single point (phi), inputs are parameters
        return np.math.factorial(p)/(np.math.factorial(i)*np.math.factorial(p-i))*(phi**i)*((1-phi)**(p-i))
    
    def C(self, phi, N1, N2):
        return phi**N1 * (1-phi)**N2
    
    def S(self,phi,P):
        sum=0
        l=len(P)-1
        for i in range(l+1):
            sum += P[i]*self.Bernstein(i,l,phi)
        return sum
    
    def zeta(self, phi, P, N1, N2):  
        return self.C(phi, N1, N2)  *  self.S(phi, P) #+  phi*.0005
    
    def CST_curve(self , P, points):
        #defines coordiantes of points which define CST curve
        
        self.x = (np.logspace(0,1, base = 10, num = points)-1)/9
        
        for i in range(points):
            #self.x.append( i/100. )
            self.y.append( self.zeta( self.x[i], self.P, self.N1, self.N2) )   
        self.y = np.asarray(self.y)
   
    def plot_CST(self, tex = False):
        #from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
        fig, ax = plt.subplots(figsize = (4,2), dpi = 300)
        ax.plot(self.x, self.y, 'k-', linewidth = 1)
        plt.grid(True)
        plt.grid(which='major', linewidth = 0.2)
        plt.title('CST curve formed using given vector P')
        plt.show()
        
        if tex:
            X = np.append(self.x.reshape(-1,1), self.y.reshape(-1,1), axis = 1)
            np.savetxt(r'E:\propeller\python\wing3d\tex-plots\cstcurve.txt', X)

    def plotBernBase(self, n):
        """ just to plot bernstein basis, for refernce """
        def Bernstein(i,p,phi):
            # returns value of bernstein for single point (phi), inputs are parameters
            return np.math.factorial(p)/(np.math.factorial(i)*np.math.factorial(p-i))*(phi**i)*((1-phi)**(p-i))
        fig, ax = plt.subplots(figsize = (4,2), dpi = 300)
        x = np.linspace(0,1,100)
        for i in range(n+1):
            print('i, p = {}, {}'.format(i,n-i))
            y =  Bernstein(i, n, x)
            ax.plot(x, y, 'k-', linewidth = 1)
            X = np.append(x.reshape(-1,1), y.reshape(-1,1), axis = 1)
            np.savetxt(r'E:\propeller\python\wing3d\tex-plots\bernstein{}.txt'.format(i), X)