# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 08:45:46 2019

@author: WK5521
N for leading edge, M for trailing edge
N1: derivative -> 0;  -> 0.5 -> infty, >0.5 shit, keep it below 1, but > 0.7
N2: want smooth LE -> N2 = 0.5, want sharp LE -> N2 = 1
M1: want derivative 0 -> M1 = 1, icreasae derivative by setting M1 -> 0, direct derivative control with M1 = 0
M2: keep at 0.5.
    


P = [0.03, .08,  .12 , .16,  .06]
Q = [.16, .18, .2,  .11, .1]

#pars = [1. 0. 1. 1. 1. 0. 0. 0. 0. 1.]
#P = pars[:5]
#Q = pars[5:]

airfoil = CSTairfoil(P, Q, N1 = .5, N2 = .5,  M1 = 1, M2 = 1, points = 150)

a = airfoil.genAirfoil(t = 0.008)

#a.runXFOIL(alfa = 2, re = 1e6)
a.plotAirfoil(saveFig = True, dpi = 200)




a.findCamberThickness(True, True, 'cst0')
a.plotAirfoil(tex = True, nametex = 'cst0')
a.runXFOIL(alfa = 2, re = 1e5, cp = True, n_crit = 2)
a.plotCp()


a.runFluent(5,.4,.5, mesh = 'o')

a1, cd, cm , cl  = a.runXFOIL(.8, re = 5e5, m = .1)
print(a1, cd, cm , cl)
a1, cd, cm , cl  = a.runXFOIL(.4, re = 1e6, m = .2)
print(a1, cd, cm , cl)



========================================================
parameters for optimzation:
max camber : P = [0.05, .3,  0.3,   .3 , .3,  .2]
min thickenss: Q = [.1, .1, .1, .1,  .1, .1]
max tthickness: Q = [.25, .25, .25, .25,  .25, .25]

"""
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("E:\propeller\python")
from airfoil import Airfoil

from subprocess import run, PIPE


class foil2(Airfoil):
    def __init__(self):
        print('created foil2')

class CSTcurve:
    """
    P = [1,1,1,1,1]
    a = CSTcurve(P, 1,1, 100)
    a.plot_CST()
    
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
        from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
        fig, ax = plt.subplots(figsize = (4,2), dpi = 300)
        ax.plot(self.x, self.y, 'k-', linewidth = 1)
        plt.grid(True)
        plt.grid(which='major', linewidth = 0.2)
        plt.title('CST curve formed using given vector P')
        plt.show()
        
        if tex:
            X = np.append(self.x.reshape(-1,1), self.y.reshape(-1,1), axis = 1)
            np.savetxt(r'E:\propeller\python\wing3d\tex-plots\cstcurve.txt', X)

def plotBernBase(n):
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
    
        
class CSTairfoil:
    """
    initializes airfoil with specified parameters
    using CST curves for top and bottom surfaces
    """
    
    def __init__(self, P, Q, N1 = 0.5 , N2 = 0.5, M1 = 1, M2 = 1, points = 160, N = None, verbose = False):
        """
        P and Q are vectors of variables for top and bottom airfoil curves
        Ns and Ms are parameters for defining smoothness of LE and TE
        points is number of points for each curve
        
        N is used for Panukl 3d wing generation, where N stands for section ID.
        
        """
        self.maxthickness = None
        self.area = None
        self.LD = None
        
        self.HFlift = None
        self.HFdrag = None
        self.LFlift = None
        self.LFdrag = None
        
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
        
        self.xfoilfile = r'E:/wing/airfoil/'  + self.name + '.txt'
        self.xfoilthickness = r'E:\propeller\XFOIL\thicknessfoil.txt' # foil with thickness distribution
        self.xfoilcamber = r'E:\propeller\XFOIL\camberfoil.txt'     # foil with camber distribution
        self.cambertxt = 'camberdist.txt'     # file with camber distribution
        self.xfoil = r'E:\propeller\XFOIL\cstfoil.txt'              # final foil
        
        self.camber = CSTcurve(P, N1, M1, points)
        self.thick = CSTcurve(Q, N2, M2, points)
        
        topy = self.camber.y + np.array(self.thick.y)/2
        boty = self.camber.y - np.array(self.thick.y)/2
        
        self.top = CSTcurve(P, N1, M1, points)
        self.bottom = CSTcurve(Q, N2, M2, points)
        self.top.y = topy
        self.bottom.y = boty
        
        
        self.addCamber_XFOIL()
    
    def plotAirfoil(self):
        
        fig = plt.figure(figsize = (6,3),dpi = 300)
        ax1   = fig.add_subplot(221)

        ax1.plot(self.top.x, self.top.y, 'k-', label = 'top', linewidth = 1)
        ax1.plot(self.bottom.x, self.bottom.y, 'k-', label = 'top', linewidth = 1)
        ax1.plot(self.camber.x, self.camber.y, 'k--', label = 'camber', linewidth = 1)
        plt.axis('equal')
#        ax1.grid()
        fig.tight_layout()
        
        
        
    def genAirfoil(self, t):
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
            file= self.xfoilfile
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
                
plotParsBounds = False
if plotParsBounds:
    
    P = [0.3, .25,  0.25,   .15 , .1,  .05]
    Q = [.3, .4, .3, .3,  .2, .2]
    airfoil = CSTairfoil(P, Q, N1 = .8, N2 = .5,  M1 = 0, M2 = .5, points = 250)
    a1 = airfoil.genAirfoil(t = 0.008)
    
    P = [0.3, .25,  0.25,   .15 , .1,  .05]
    Q = [.1, .08, .08, .05,  .05, .03]
    airfoil = CSTairfoil(P, Q, N1 = .8, N2 = .5,  M1 = 0, M2 = .5, points = 250)
    a0 = airfoil.genAirfoil(t = 0.008)
    a0.findCamberThickness(True)
    a1.findCamberThickness(True)
    
    plt.figure(figsize = (6,3),dpi = 200)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.plot(a0.camber[:,0], a0.camber[:,1], 'k-',linewidth = 1.2, label = 'max camber')
    plt.plot( a0.thickness[:,0], a0.thickness[:,1], 'k--',linewidth = 1, label = 'min thickness')
    plt.plot( a0.thickness[:,0], -a0.thickness[:,1], 'k--',linewidth = 1)
    plt.plot( a1.thickness[:,0], -a1.thickness[:,1], 'k-.',linewidth = 1, label = 'max thickness')
    plt.plot( a1.thickness[:,0], a1.thickness[:,1], 'k-.',linewidth = 1)
    
    plt.xlabel(r'$x/c$',fontsize=12)
    plt.ylabel(r'$y/c$',fontsize=12)
    plt.title(r"Optimization bounds", fontsize=12)
    # Make room for the ridiculously large title.
    plt.subplots_adjust(top=0.8)
    plt.axis('equal')
    plt.legend()
    plt.tight_layout()
#        plt.grid(True)
    saveFig = True
    if saveFig:
        plt.savefig(r'C:\Users\wk5521\Documents\python\saved_plots\bounds', dpi = 1000)
    plt.show()
    
    
    
class RunDesign:
    """ 
    
    design = RunDesign(r'H:\airfoil design\omesh/')
    design.runSample()
    design.runGP()
    
    
    design.runKriging()
    
    """
    def __init__(self, path):
        
        pars = 11
        samples = 100
        
        self.path = path # main path to folder
        self.X = np.loadtxt(r'{}X,{},{}.txt'.format(path, pars, samples))
        self.bounds = np.loadtxt(r'{}bounds.txt'.format(path))
        self.Ubnds = self.bounds[:,1]
        self.Lbnds = self.bounds[:,0]
        self.objective = []
        
        # defines settings
        self.mach = .5
        self.chord = .5
        self.Cl = np.zeros((np.shape(self.X)[0],))
        self.Cd = np.zeros((np.shape(self.X)[0],))
        self.Cm = np.zeros((np.shape(self.X)[0],))
        self.t = np.zeros((np.shape(self.X)[0],))
        
        
    def runSample(self):
#        self.Cl = np.zeros(np.shape(self.X))
#        self.Cd = np.zeros(np.shape(self.X))
#        self.Cm = np.zeros(np.shape(self.X))
       
        try:
            self.t = np.loadtxt(r'{}t.txt'.format(self.path))
            
        except IOError:
            i=0
            for case in self.X:
                self.t[i] = self.runPointThickness(case, [.3])
                i+=1
            np.savetxt(r'{}t.txt'.format(self.path), self.y)
        i=0
        for case in self.X:
            self.Cl[i], self.Cd[i] , self.Cm[i] = self.runPoint(case, self.mach)
#            self.t[i] = self.runPointThickness(case, [.3])
            i+=1
            
            
    def returnAirfoil(self, x):
        num = int((len(x)-1)/2)
        P = self.Lbnds[:num] + self.Ubnds[:num] * x[:num]
        Q = self.Lbnds[num:-1] + self.Ubnds[num:-1] * x[num:-1]
        Q[0] = P[0] + self.Lbnds[num] + self.Ubnds[num]*x[num]
        
        airfoil = CSTairfoil(P, Q, N1 = .5, N2 = .5,  M1 = 1, M2 = 1, points = 150)
        
        a = airfoil.genAirfoil(t = 0.01)
        return a


    def runPoint(self, x, mach, rho=1.225, visc=1.78e-5):
        """  method to run point calculation for defined x 
        x[-1] is aoa  """

        aoa = self.Lbnds[-1] + self.Ubnds[-1] * x[-1]
        
        name = ''
        for num in x[:-1]:
            name += '{:.3f},'.format(num)
        name += '{:.2f},'.format(aoa)
        name += '{:.2f}'.format(mach)
        
        try:
            c = np.loadtxt('{}reports/{}.out'.format(self.path, name), skiprows = 200)
            cl = (c[-1, 1]+c[-2, 1])/2
            cd = (c[-1, 2]+c[-2, 2])/2
            cm = (c[-1, 3]+c[-2, 3])/2
            print('loaded data from file')
            return cl, cd, cm
        
        except IOError:
            
            num = int((len(x)-1)/2)

            x = np.asarray(x)
            
            P = self.Lbnds[:num] + self.Ubnds[:num] * x[:num]
            Q = self.Lbnds[num:-1] + self.Ubnds[num:-1] * x[num:-1]
            Q[0] = P[0] + self.Lbnds[num] + self.Ubnds[num]*x[num]
            
            airfoil = CSTairfoil(P, Q, N1 = .5, N2 = .5,  M1 = 1, M2 = 1, points = 150)
            
            a = airfoil.genAirfoil(t = 0.01)
            a.plotAirfoil(saveFig = True, dpi = 200)
            
            a.runFluent(aoa, mach, self.chord, name = name, path = self.path, rho = rho, viscosity = visc, mesh = 'o')
            
            reportpath = self.path + 'reports/{}.out'.format(name)
            
            c = np.loadtxt(reportpath, skiprows = 200)
            cl = (c[-1, 1]+c[-2, 1])/2
            cd = (c[-1, 2]+c[-2, 2])/2
            cm = (c[-1, 3]+c[-2, 3])/2
#            print('cl, cd, cl/cd:  {:.3f},{:.3f},{:.1f}'.format(cl, cd, cl/cd))
            return cl, cd, cm
        

    def runPointThickness(self, x, x_c):
        """
        inputs are x and [x_c] which is a position(s) of required thickness(es)
        
        currently supports single thickness only
        
        """

        num = int((len(x)-1)/2)
        x = np.asarray(x)
        
        P = self.Lbnds[:num] + self.Ubnds[:num] * x[:num]
        Q = self.Lbnds[num:-1] + self.Ubnds[num:-1] * x[num:-1]
        Q[0] = P[0] + self.Lbnds[num] + self.Ubnds[num]*x[num]
        
        airfoil = CSTairfoil(P, Q, N1 = .5, N2 = .5,  M1 = 1, M2 = 1, points = 150)
        
        a = airfoil.genAirfoil(t = 0.01)
        
        for t in x_c:
            t = a.t_x(t)
        return t
    
    def saveGPScales(self, path = None):
        if path is None:
            path = r'C:\Users\wk5521\OneDrive\doktorat\beta\results\airfoil\length_scales.txt'
            
        a1 = 1/self.gp_Cl.kernel_.get_params()['k2__length_scale'].reshape(-1,1)
        a2 = 1/self.gp_Cd.kernel_.get_params()['k2__length_scale'].reshape(-1,1)
        a3 = 1/self.gp_Cm.kernel_.get_params()['k2__length_scale'].reshape(-1,1)
        a4 = 1/self.gp_t.kernel_.get_params()['k2__length_scale'].reshape(-1,1)
        A = np.append(a1, a2)
        A = np.append(a1, a2, axis = 1)
        A = np.append(A, a3, axis = 1)
        A = np.append(A, a4, axis = 1)
        np.savetxt(path, A)
        
    def runGP(self, plot=False, plot3d = False, dpi = 200):
        
        
        from sklearn.gaussian_process import GaussianProcessRegressor, GaussianProcess
        from sklearn.gaussian_process.kernels import RBF, Matern, ConstantKernel as C
        from scipy.optimize import differential_evolution as de    
        
#        def minim(x, cl, cd):
#            a = self.returnAirfoil(x)
#            t = a.t_x(.3) 
#            pen = .12 - t if t < .12 else 0
#            return -cl / (cd + pen)

#        y = np.loadtxt('H:/airfoil design/omesh/y.txt')
        
        kernel = C(1.0, (1e-3, 1e3)) * Matern(length_scale=[1.0]*11, length_scale_bounds=(1e-1, 100.0),nu=5/2) # think about appropriate kernel
        
        gp_Cl = GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer=9)
        gp_Cd = GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer=9)
        gp_Cm = GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer=9)
        gp_t = GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer=9)
        gp_Cl.fit(self.X, self.Cl)
        gp_Cd.fit(self.X, self.Cd)
        gp_Cm.fit(self.X, self.Cm)
        gp_t.fit(self.X, self.t)
        self.gp_Cl = gp_Cl
        self.gp_Cd = gp_Cd
        self.gp_Cm = gp_Cm
        self.gp_t = gp_t    
        
        def f(x):
            # objective function maximizing l/d
            if self.gp_t.predict(x.reshape(1,-1)) < .12:
                return -self.gp_Cl.predict(x.reshape(1,-1)) / (self.gp_Cd.predict(x.reshape(1,-1))+self.gp_t.predict(x.reshape(1,-1)))
            else:
                return -self.gp_Cl.predict(x.reshape(1,-1)) / (self.gp_Cd.predict(x.reshape(1,-1)))
            
        def f2(x):
            # objective function minimizing cd with cl constraint-penalty
            if self.gp_Cl.predict(x.reshape(1,-1)) < 1 and self.gp_t.predict(x.reshape(1,-1)) < .12:
                return (self.gp_Cd.predict(x.reshape(1,-1)) + (1-self.gp_Cl.predict(x.reshape(1,-1))) ) + 1/self.gp_t.predict(x.reshape(1,-1))
            elif  self.gp_Cl.predict(x.reshape(1,-1)) < 1: 
                return self.gp_Cd.predict(x.reshape(1,-1)) + (1-self.gp_Cl.predict(x.reshape(1,-1))) 
            elif self.gp_t.predict(x.reshape(1,-1)) < .12:
                return self.gp_Cd.predict(x.reshape(1,-1)) + 1/self.gp_t.predict(x.reshape(1,-1))
            else:
                return self.gp_Cd.predict(x.reshape(1,-1))
            
        f = f2
        
        for case in self.X:
            self.objective.append(f(case))
        
        bounds = [(0,1)]*11
        result = de(f , bounds = bounds)
        print(result)
        
        cl, cd, cm = self.runPoint(result.x, self.mach)
        self.X = np.append(self.X, result.x.reshape(1,-1), axis = 0)
        self.Cl = np.append(self.Cl, cl)
        self.Cd = np.append(self.Cd, cd)
        self.Cm = np.append(self.Cm, cm)
        self.t = np.append(self.t, self.runPointThickness(result.x, [.3]))
        
#        print(gp_Cl.kernel_)
#        print(gp_Cd.kernel_)
#        print(gp_Cm.kernel_)
        
        for i in range(5):
            print('\n\ninfill {}'.format(i))
#            print(result.x)
#            print(self.X.shape)
#            y = np.append(y, -cl/cd)
            
            print('fitting GPs')
            gp_Cl.fit(self.X, self.Cl)
            gp_Cd.fit(self.X, self.Cd)
            gp_Cm.fit(self.X, self.Cm)
            gp_t.fit(self.X, self.t)
            
#            gp = GaussianProcessRegressor(kernel = kernel, n_restarts_optimizer=9)
#            gp.fit(self.X, y)
#            def f(x):
#                return gp.predict(x.reshape(1,-1))
            result = de(f , bounds = bounds)
            
            clp = self.gp_Cl.predict(result.x.reshape(1,-1))
            cdp = self.gp_Cd.predict(result.x.reshape(1,-1))

            
            
            cl, cd, cm = self.runPoint(result.x, self.mach)
            print('predicted cl, cd, cl/cd: {:.3f} {:.3f} {:.1f}'.format(clp[0], cdp[0], clp[0]/cdp[0]))
            print('actual cl, cd, cl/cd: {:.3f} {:.3f} {:.1f}'.format(cl, cd, cl/cd))
            self.objective.append(f(result.x))
            
            self.X = np.append(self.X, result.x.reshape(1,-1), axis = 0)
            self.Cl = np.append(self.Cl, cl)
            self.Cd = np.append(self.Cd, cd)
            self.Cm = np.append(self.Cm, cm)
            self.t = np.append(self.t, self.runPointThickness(result.x, [.3]))
        
        i = np.argmin(self.objective)        
        self.finalAirfoil = self.returnAirfoil(self.gp_Cl.X_train_[i])
        
        if plot:
            # plotting and postprocessing part
            plt.figure(dpi = 150)
            plt.plot(1/self.gp_Cl.kernel_.get_params()['k2__length_scale'], 'ks-',mfc = 'none', label = '$C_L$ GP')
            plt.plot(1/self.gp_Cd.kernel_.get_params()['k2__length_scale'], 'ko--',mfc = 'none', label = '$C_D$ GP')
            plt.plot(1/self.gp_Cm.kernel_.get_params()['k2__length_scale'], 'kh:',mfc = 'none', label = '$C_M$ GP')
            plt.plot(1/self.gp_t.kernel_.get_params()['k2__length_scale'], 'kX-.',mfc = 'none', label = '$t$ GP')
            plt.legend()
            plt.grid('major', linewidth = .1)
            plt.ylabel('inverted length scale')
            plt.xlabel ('parameter number')
            plt.tight_layout()
            plt.show()
        
        if plot3d:
            # 3d plotting
            # x1 and x2 are indices to be plotted
            
            gp = self.gp_Cl
            
            x1 = 0
            x2 = 10
            
            x_vals = np.arange(0, 1, .01)
            y_vals = np.arange(0, 1, .01)
            X, Y = np.meshgrid(x_vals, y_vals)
            
            Z = np.zeros((100, 100))
            sigma = np.zeros((100, 100))
            for i in range(100):
                for j in range(100):
                    x = np.ones(11)*0.5
                    x[x1] = x_vals[i]
                    x[x2] = y_vals[j]
                    Z[j,i], sigma[j,i] = gp.predict(x.reshape(1,-1), return_std = True)
            
            
            fig = plt.figure(dpi = dpi)
            cp = plt.contourf(X, Y, Z, cmap="Greys")
            plt.colorbar(cp)
    #        plt.clabel(cp, inline=True,  fontsize=10)
            plt.title('objective function prediction')
            plt.xlabel('$x[{}]$'.format(x1))
            plt.ylabel('$x[{}]$'.format(x2))
            plt.tight_layout()
            plt.show()
                 
            fig = plt.figure(dpi = dpi)
            cp = plt.contourf(X, Y, sigma, cmap="Greys")
            plt.colorbar(cp)
    #        plt.clabel(cp, inline=True,  fontsize=10)
            plt.title('standard error')
            plt.xlabel('$x[{}]$'.format(x1))
            plt.ylabel('$x[{}]$'.format(x2))
            plt.tight_layout()
            plt.show()
            
    def runKriging(self):
        from pyKriging.krige import kriging  
#        kl = kriging(self.X, self.Cl)
#        kl.train()
#        kd = kriging(self.X, self.Cd)
#        kd.train()
        
    
        
        
        # construct minimization objective with thickness constraints
        def minim(x, cl, cd):
            a = self.returnAirfoil(x)
            t = a.t_x(.3) 
            pen = .12 - t if t < .12 else 0
            return -cl / (cd + pen)
        
#        y = np.zeros(len(self.X))
#        for i in range(len(self.X)):
#            y[i] = minim(self.X[i], self.Cl[i], self.Cd[i])
#        np.savetxt('H:/airfoil design/omesh/y.txt', y)
        
        y = np.loadtxt('H:/airfoil design/omesh/y.txt')
        
        print('building initial kriging')
        self.k = kriging(self.X, y)
        self.k.train()
        from scipy.optimize import differential_evolution as de    
        bounds = [(0,1)]*11
        result = de(self.k.predict , bounds = bounds)
        print(result)
        
        
        for i in range(15):
            
            cl, cd, cm = self.runPoint(result.x, self.mach)
            yi = minim(result.x, cl, cd)
            self.k.addPoint(result.x, yi)
            self.k.train()
            result = de(self.k.predict , bounds = bounds)
        
        i = np.argmin(self.k.y)
        self.finalAirfoil = self.returnAirfoil(self.k.X[i])
        
        
#        def minimT(x):
#            a = self.returnAirfoil(x)
#            t = a.t_x(.3) 
#            if t < .12:
#                return -kl.predict(x) /  (kd.predict(x) + t)
#            else:
#                return -kl.predict(x) /  kd.predict(x)
#            
#            
#        
#        from scipy.optimize import differential_evolution as de    
#        bounds = [(0,1)]*11
##        self.result = de(minim, bounds = bounds)
#        
##        cl, cd, cm = self.runPoint(self.result.x, self.mach)
##        
##        print('error = {:.2f}'.format( (kl.predict(self.result.x)/kd.predict(self.result.x) - cl/cd)/(cl/cd)))
#        
##        for i in range(numiter):  
##        print 'Infill iteration {0} of {1}....'.format(i + 1, numiter)
#        newpoints = kd.infill(10)
#        for point in newpoints:
#            cl, cd, cm = self.runPoint(point, self.mach)
#            kl.addPoint(point, cl)
#            kd.addPoint(point, cd)
#        kl.train()
#        kd.train()
#
#        
#        for i in range(5):  
#            result = de(minimT, bounds = bounds)
#            clp = kl.predict(result.x)
#            cdp = kd.predict(result.x)
#            
#            cl, cd, cm = self.runPoint(result.x, self.mach)
#            print('predicted cl, cd,cl/cd= {:.3f},{:.3f},{:.1f} at aoa = {:.1f}'.format(clp, cdp, clp/cdp, result.x[-1]))
#            print('actual, cl, cd, cl/cd = {:.3f},{:.3f},{:.1f} at aoa = {:.1f}'.format(cl, cd, cl/cd, result.x[-1]))
#            kl.addPoint(result.x, cl)
#            kd.addPoint(result.x, cd)
#            kl.train()
#            kd.train()
#            print('shape of X: {}'.format(kl.X.shape))
#            
#    



