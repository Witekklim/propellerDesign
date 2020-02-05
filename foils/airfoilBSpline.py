# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 14:30:14 2019

@author: WK5521

knotsU = np.array([[0,0], [0., .04], [0.2, .1], [.4, .12],[.6,.1], [.8, .05] ,[1,0]   ])
knotsL = np.array([[0,0], [0., -.025], [0.2, -.05], [.4, -0.06],[.6, -.01],[.8, .01], [1,0]   ])
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
sys.path.append("E:\propeller\python")
from airfoil import Airfoil



class BSplineAirfoil(object):
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
    def __init__(self, knots, points = 10):
        """ finds bspline curve for specified knots
        inputs: knots coords as numpy array  [  [x0,y0],[x1,y1],...,[xn-1,yn-1]  ]
        attribute: coords of points on spline, number defined in __init__  [ [x0,x1,x2,...,xn-1], [y0,y1,y2,...,yn-1] ]
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
        
        
     
        
  





      
        
from pyKriging.samplingplan import samplingplan
def runSample(Ubnds, Lbnds, k=11, n=200, La = 1, Ha = 5, start = 0, num = 0, folder = '0'):
    
    # create sampling plan
    if False:
        sp = samplingplan(k)  
        X = sp.optimallhc(n)
        np.savetxt(r'H:\airfoil design\X,{},{}.txt'.format(k,n), X)
    
    X = np.loadtxt(r'H:\airfoil design\X,{},{}.txt'.format(k,n))
    
    if num is 0:
        num = n

    # case 117
    cl = np.zeros(n)
    cd = np.zeros(n)
    cm = np.zeros(n)
    i = 0
    for case in X[start:start+num]:
#        runPoint(case, Ubnds, Lbnds, La, Ha, folder)
        cl[i], cd[i], cm[i] = runPoint(case, Ubnds, Lbnds, La, Ha, folder)
        i+=1


def runPoint(x, Ubnds, Lbnds, La , Ha, folder = '0', explicit = True):
    """ La , Ha are lower and upper aoa bounds  """
    name = ''
    for num in x:
        name += '{:.3f},'.format(num)
    
    u = Ubnds[:,0] + x[:5] * (Ubnds[:,1] - Ubnds[:,0])
    l = Lbnds[:,0] + x[5:10] * (Lbnds[:,1] - Lbnds[:,0])
    aoa = La + Ha* x[10]
    name += '{:.3f}'.format(aoa)
    
    print(name)

    
    try:
        c = np.loadtxt('H:/airfoil design/reports{}/{}-alfa{:.2f}-mach-0.20.out'.format(folder, name,aoa ), skiprows = 500)
        cl = (c[-1, 1]+c[-2, 1])/2
        cd = (c[-1, 2]+c[-2, 2])/2
        cm = (c[-1, 3]+c[-2, 3])/2
        print('returning data from output')
        return cl, cd, cm
    except OSError:
    
        knotsU = np.array([[0,0], [0., u[0]], [0.2, u[1]], [.4, u[2]],[.6, u[3]], [.8, u[4]],[1,0]   ])
        
        if explicit:
            knotsL = np.array([[0,0], [0., l[0]], [0.2, l[1]], [.4, l[2]],[.6, l[3]], [.8, l[4]],[1,0]   ])
        else:
            knotsL = np.array([[0,0], [0., u[0] - l[0]], [0.2, u[1] - l[1]], [.4, u[2] - l[2]],[.6, u[3] - l[3]], [.8, u[4] - l[4]],[1,0]   ])
    
        foil = BSplineAirfoil(knotsU, knotsL)
        foil.plotAirfoil(True)
        a = foil.genAirfoil()
        
        a.runFluent(aoa, .2, 1, name, folder)
        
        c = np.loadtxt('H:/airfoil design/reports{}/{}-alfa{:.2f}-mach-0.20.out'.format(folder, name, aoa), skiprows = 500)
        cl = (c[-1, 1]+c[-2, 1])/2
        cd = (c[-1, 2]+c[-2, 2])/2
        cm = (c[-1, 3]+c[-2, 3])/2
        return cl, cd, cm
    
    
def returnT(x, Ubnds, Lbnds):
    u = Ubnds[:,0] + x[:5] * (Ubnds[:,1] - Ubnds[:,0])
    l = Lbnds[:,0] + x[5:10] * (Lbnds[:,1] - Lbnds[:,0])
    
    knotsU = np.array([[0,0], [0., u[0]], [0.2, u[1]], [.4, u[2]],[.6, u[3]], [.8, u[4]],[1,0]   ])
    knotsL = np.array([[0,0], [0., l[0]], [0.2, l[1]], [.4, l[2]],[.6, l[3]], [.8, l[4]],[1,0]   ])
    
    foil = BSplineAirfoil(knotsU, knotsL)
    foil.plotAirfoil(True)
    a = foil.genAirfoil()
    return a.t_x(.3)
        
def createKriging(Ubnds, Lbnds, La , Ha , folder = '3', reqcl = 0):
    from pyKriging.krige import kriging  

    # collect data
    
    X = np.loadtxt(r'H:\airfoil design\X,11,100.txt')
    t = np.loadtxt(r'H:\airfoil design\t,11,100.txt')
#    Ubnds = np.array([ [0.02, .04], [0.07, .1],   [.065, .1],  [.05,.09],   [.01, .06]   ])
#    Lbnds = np.array([ [-.01, -.03], [-.025, -.05], [-.02, 0.0],[-.02, .01], [-.005, .01] ])
    num = len(X)
    
    
    cl = np.zeros(num)
    cd = np.zeros(num)
#    t = np.zeros(num)
    i = 0
    for case in X:
        print('collecting data for case {}'.format(i))
        name = 'H:/airfoil design/reports{}/'.format(folder)
        for num in case:
            name += '{:.3f},'.format(num)
        aoa = 1 + 6 * case[10]
        name += '{:.3f}-alfa{:.2f}-mach-0.20.out'.format(aoa,aoa)
        c = np.loadtxt(name, skiprows = 500)
        cl[i] = (c[-1, 1]+c[-2, 1])/2
        cd[i] = (c[-1, 2]+c[-2, 2])/2
#        t[i] = returnT(case, Ubnds, Lbnds)
        i += 1
    
    if False:
        np.savetxt(r'H:\airfoil design\t,11,100.txt', t)
    
    # maximization of cl/cd
    
    cls = np.asarray(cl)
    cds = np.asarray(cd)
    ts = np.asarray(t)
    
    def merit(reqcl, cl, cd, ts):
        if reqcl == 0:
            return  cl/cd
        else:
            # return penalized drag by scaled value of cl constraint breach
            aero =  np.where(cl >=reqcl, cd, cd+(reqcl - cl))
            total = np.where(ts > .11, aero, aero + (.11- ts))
            return total

    y = merit(reqcl, cls, cds, ts)
    
    # Now that we have our initial data, we can create an instance of a Kriging model
#    k = kriging(X, y, name='simple')  
#    k.train()
    
    # create two kriging models: for cl and cd separately
    kl = kriging(X, cls)
    kl.train()
    kd = kriging(X, cds)
    kd.train()
    kt = kriging(X, ts)
    kt.train()
    
    # this part is used to improve surrogate global quality
    numiter = 0 
    for i in range(numiter):  
        print ('Infill iteration {0} of {1}....'.format(i + 1, numiter))
        newpoints = k.infill(5)
        for point in newpoints:
            print('exploration')
            cl , cd, cm = runPoint(point, Ubnds, Lbnds, La, Ha , folder = folder)
            k.addPoint(point, merit(reqcl, cl, cd))
        k.train()
    
    from scipy.optimize import differential_evolution as de
    
    n = 3
    for i in range(n):
        print('exploitation {} of {}'.format(i+1, n))
        
        def minim(x, *args):
            reqcl = args[0]
            cl = kl.predict(x)
            cd = kd.predict(x)
            t = kt.predict(x)
            re = merit(reqcl, cl, cd, t)
            if cl>reqcl:
                res = cd
            else:
                res = cd + (reqcl-cl)
            return re
        
        args = ([reqcl])
        bounds = [(0,1)]*11
        result = de(minim, bounds = bounds, args = args )
        print(result)
        
        cl , cd, cm = runPoint(result.x, Ubnds, Lbnds, La, Ha , folder = folder)
        t = returnT(result.x, Ubnds, Lbnds)
        print('calculated values cl: {}, cd: {} cl/cd: {}'.format(cl, cd, cl/cd))
        
        kl.addPoint(result.x, cl)
        kl.train()
        kd.addPoint(result.x, cd)
        kd.train()
        kt.addPoint(result.x, cd)
        kt.train()
    
  
 
#createKriging(Ubnds, Lbnds, La, Ha)
     
run = False
if run:
    Ubnds = np.array([ [0.02, .06], [0.07, .12],   [.065, 0.12],  [.05,.08],   [.01, .05]   ])
    Lbnds = np.array([ [-.015, -.03], [-.02, -.06], [-.015, -.025],[.005, -.005], [.005, -.02,] ])
            
    La = 1
    Ha = 6   
        
    #runSample(Ubnds, Lbnds, k=11, n=100, La = La , Ha= Ha , start = 0, num = 0, folder = '3')
        
    createKriging(Ubnds, Lbnds, La, Ha, folder = '3', reqcl = .6)
        
    # new      
    
    La = 1
    Ha = 6   
        
    
    Ubnds = np.array([ [0.02, .06], [0.07, .12],   [.065, 0.12],  [.05,.08],   [.01, .05]   ])
    Lbnds = np.array([ [-.015, -.03], [-.02, -.06], [-.015, -.02],[.005, -.005], [.005, -.02,] ])
    
    
#    runSample(Ubnds, Lbnds, k=11, n=100, La=La, Ha=Ha, start = 2, num = 0, folder = '4')
    
    createKriging(Ubnds, Lbnds, La, Ha, folder = '4', reqcl = .6)
            
            

















