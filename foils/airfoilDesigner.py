# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 08:50:22 2019

@author: Witold Klimczyk
"""

from airfoilCST import CSTairfoil
from airfoilBSpline import BSplineAirfoil
from pyKriging.samplingplan import samplingplan
from scipy.optimize import minimize, differential_evolution as de
import numpy as np

class airfoilDesigner():
    """  used to:
            kriging based design
            differential evolution based design   
    """
    def __init__(self, L_req, parametrization , N_pars = None,  constraints = None ,
                 ts = None, xs = None, weights = None):
        """
        L_req: array of cls for multi-point design: e.g. [.3,.6,.9]
        parametrization: 'CST' or 'BS' (b-spline)
        """
        print('running airfoil designer')

        self.parametrization = parametrization
        self.L_req = L_req
        self.N_pars = N_pars
        
        # by default top and bottom/camber and thickness curves are of the same size
        self.N_1 = int(N_pars/2)
        self.N_2 = int(N_pars/2)
        
        self.cm_max = []
        
        self.tBounds = None
        self.cBounds = None
        
        self.fileSampling = r'E:\propeller\python\wing3d\X-{}-{}.txt'.format(N_pars, self.N_samples)
        self.filey = r'E:\propeller\python\wing3d\y-{}-{}.txt'.format(N_pars, self.N_samples)
        
        # constraints
        if ts is None:
            self.ts = np.array([.11, .08, .03])
            self.xs = np.array([.3, .5, .8])
            self.weights = np.array([1.,1.,1.]) # defines weightning factor for breached constraints 
        else:
            self.ts = ts
            self.xs = xs
            self.weights = weights
        # observations
        self.y = np.zeros(self.N_samples)
        
        # bounds section
        self.LB= None
        self.UB = None
        self.bounds = None
        self.init_bounds()
        
        self.boundsBS = np.array([ [0.02 , 0.040] , [0.06, 0.1] , [0.06, 0.10] , [0.05, 0.1] , [0.02, 0.05]
                                  ,[0.045 ,0.065] , [0.09,0.16] , [0.08, 0.14] , [0.06, 0.1] , [0.03, 0.08] ])
        
        self.finalAirfoilBS = None
        self.finalAirfoil = None
    
    def airfoilfromX(self, X):
        X_real = X * (self.UB - self.LB) + self.LB
        airfoil = CSTairfoil( X_real[:self.N_1], X_real[self.N_1:], N1 = .5, N2 = .5,  M1 = 1, M2 = 1 )
        airfoil.addCamber_XFOIL()
        a = airfoil.genAirfoil(t = 0.003)
        return a
    
    
    def airfoilBSfromX(self, x):
        x1 = np.zeros(len(x))
        for i in range(len(x)):
            x1[i] = self.boundsBS[i,0] + x[i] * ( self.boundsBS[i,1] - self.boundsBS[i,0] ) 
        x = x1
        knotsU = np.array([[0,0], [0., x[0]],         [0.2, x[1]],        [.4, x[2]],       [.6, x[3]],        [.8, x[4]],        [1,0]   ])
        knotsL = np.array([[0,0], [0., x[0] - x[5] ], [0.2, x[1] - x[6]], [.4, x[2] - x[7]],[.6, x[3] - x[8]], [.8, x[4] - x[9]], [1,0]  ])
        foil = BSplineAirfoil(knotsU, knotsL)
        a = foil.genAirfoil(t = 0.003)
        return a
    
    
    def init_bounds(self):
        """ initialize bounds for CST based optimization   """
        # camber
        print(  'initializing bounds')

        c_LB = np.array(np.ones(self.N_2)) * 0.
        c_UB = np.array(np.ones(self.N_2)) * .2
        c_UB[ 0 ] = .05

        # thickness
        t_LB = np.array(np.ones(self.N_2))*.05
        t_UB = np.array(np.ones(self.N_2))*.25
        t_LB[0] = .08
        
        LB = np.concatenate((c_LB, t_LB))
        UB = np.concatenate((c_UB, t_UB))
        
        bounds = []
        for i in range(len(UB)):
            bounds.append( (LB[i] , UB[i] ) )
        self.bounds = bounds
        self.LB = LB
        self.UB = UB
        
        
    def runOptimization(self, maxiter = 10):
#        self.result = minimize(self.merit, x0 = self.UB, method='L-BFGS-B', bounds = self.bounds, options = {'eps' : 1e-2})
        
        self.result = de(self.merit, bounds = [(0,1)] * self.N_pars, maxiter=maxiter)
        
        if self.parametrization == 'CST':
            self.finalAirfoil = self.airfoilfromX(self.result.x)
        else:
            self.finalAirfoil = self.airfoilBSfromX(self.result.x)
        self.finalAirfoil.plotAirfoil()
        
    
    def merit(self, X):
        if self.parametrization == 'CST':
            a = self.airfoilfromX(X)
        else:
            a = self.airfoilBSfromX(X)
        
        a.findCamberThickness()
        # obtain BB model for created airfoil 
        merit = 0
        for cl in self.L_req:
            cl, re, m = (cl, 1e6, .2)
            alfa, cd, cm , cl = a.runXFOIL(cl = cl , re = re, m = m, n_crit = 6)
            if (alfa, cd, cm, cl) == (1,1,1,1):
                merit = .03
            else:
                merit += cd
        
        # apply penalties for constrained designs
        for i in range(len(self.cm_max)):
            if cm > self.cm_max[i]:
                merit += .1* (cm - self.cm_max)
        for i in range(len(self.ts)):
            try:
                tx = a.t_x(self.xs[i])
                if tx < self.ts[i]:
                    merit += ( self.ts[i] - tx ) * self.weights[i]
            except TypeError:
                print('invalid thickness')
                merit = .03
        return merit

        

