# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:23:10 2019

@author: WK5521

script to run qprop automatically for different propellers, eventually use it for finding optimal design


#===================================================
the general way of of defining rotor


betas = np.linspace(20, 5, 15)
chords = np.linspace(.2, .05, 15)
rs = np.linspace(0.2, 0.95, 15)
ts = np.linspace(0.025,.005,15)

rotor = Rotor(rpm = 1500, betas = betas, rs = rs, chords = chords, ts = ts)



"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from airfoil import Airfoil

from scipy.optimize import minimize, differential_evolution, fmin
from scipy import interpolate

           
            
def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return proc_stdout
        

class Rotor:
    def __init__(self, rs = None, chords = None, betas = None, rpm = 2000, ts = None,
                 show = True, saveOutput = True, name = None, V_in = 0, N = 2,
                 ):
        """inputs:N
                D(iameter) meters
                rpm
                betas: numpy array of betas
                chord: numpy array of chords
                r: numpy array of radial positions
                t: thickness constraints (absolute)
    
        """
   
        self.D = rs[-1]*2  # diameter 
        self.rpm = rpm
        
        self.V_in = V_in # forward velocity
        self.N_stations = len(chords)
        self.rs = rs
        
        # number of blades
        self.N = N 
        
        # stations along radius
        self.chords = chords
        self.betas = betas

        # thickness constraint
        try:
            self.ts = ts / self.chords
        except TypeError:
            self.ts = None

        if show:
            self.plots()
            
        self.filename = 'E:/propeller/qprop/cases/input'
        self.motorname = 'E:/propeller/qprop/cases/motor'
        self.foilProps()
        self.CreateDataFile2()
        self.CreateMotorFile()
        self.runAnalysis(saveOutput, name)
        
        if show:
            print(self.results)
            
            
    def foilProps(self):
         # -------------------------------------------------------------------------------
        # section to define airfoils and local properties
        # -------------------------------------------------------------------------------
        def foilname(t):
            """ return airfoil name for given thickness, in % chord """
#            thicknesses = np.array([9.8, 11, 13, 14.7, 16.2])
#            names = ['mh117','mh115', 'mh114','mh113','mh112' ]
            
            thicknesses = np.array([.05])
            names = ['mh117']
            return 'mh117'
            
            a = np.insert(thicknesses, 0, 0)
            thresholds = thicknesses - 0.5 * (thicknesses - a[:-1])
            for i in range(len(thresholds)-1):
                if t < thresholds[i+1] and t>thresholds[i]:
                    return names[i]
            return names[i+1]
        
        self.airfoils = []
        self.airfoilsdata = [] # use this to load preevaluated airofil data
        # list of airfoil objects for each radial postion
        
        for i in range(self.N_stations):
#            print(foilname(self.t[i] / self.chord[i]*100))
            self.airfoils.append(Airfoil('XFOIL', r'E:/propeller/mh_airofils/{}/{}xfoil.txt'.format(foilname(self.ts[i] *100), foilname(self.ts[i] * 100) ), origin = 0))
        
        for name in ['mh112','mh113','mh114','mh115','mh117']:
            self.airfoilsdata.append(np.loadtxt(r'E:\propeller\qprop\airfoildata\{}.txt'.format(name)))
        def returnX(name):
            i = 0
            for nm in ['mh112','mh113','mh114','mh115','mh117']:
                if nm == name:
                    return self.airfoilsdata[i]
                i += 1
        
        # arrays of airofil properties
        self.Cl0 = np.zeros(self.N_stations)
        self.Cla = np.zeros(self.N_stations)
        self.Cd0 = np.zeros(self.N_stations)
        self.ClCd0 = np.zeros(self.N_stations)
        self.Cd2u = np.zeros(self.N_stations)
        self.Cd2l = np.zeros(self.N_stations)
        
        self.m = self.rs * 2 * self.rpm * np.pi / 60. / 340.
        self.re = self.m * 340 * self.chords / 1.42e-5
        
        print('\nr [m]\tchord\tfoil\tmach\tRe')
        for i in range(self.N_stations):
            print('{:.2f}\t{:.2f}\t{}\t{:.2f}\t{:.1f}'.format(self.rs[i],self.chords[i], foilname(self.ts[i] *100 ), self.m[i], self.re[i]  ) )
        
        for i in range(self.N_stations):
            # first option when dont have airfoil data
            # second uses preevaluated airofil polar data
            
#            self.Cl0[i], self.Cla[i], self.Cd0[i], self.ClCd0[i], self.Cd2u[i], self.Cd2l[i] = self.airfoils[i].qpropData(self.m[i], self.re[i])
            self.Cl0[i], self.Cla[i], self.Cd0[i], self.ClCd0[i], self.Cd2u[i], self.Cd2l[i] = readQPolar(self.re[i], X = returnX(foilname(self.ts[i] / self.chords[i]*100 )))
        
        
        self.Clmin = np.ones(self.N_stations)*.2
        self.Clmax = np.ones(self.N_stations)* 1
        self.Reref = self.chords * self.rs

    def returnMerit(self):
        """ returns merit function   """
        
        if self.Q > 70.:
            return self.P + 10000*(self.Q-70)
        if self.T < 630.:
            return self.P + 1000*(630-self.T)
                   
        return self.P
        
    def runAnalysis(self, saveOutput =False, name = None):
#        cmd = 'E:/propeller/qprop/bin/qprop '+ self.filename + ' ' + self.motorname
        if saveOutput:
            if name is None:
                subprocess_cmd('E:/propeller/qprop/bin/qprop E:/propeller/qprop/cases/input E:/propeller/qprop/cases/motor {} {} > E:/propeller/qprop/cases/output.txt'.format(self.V_in, self.rpm))
            else:
                print('saving under name {}'.format(name))
                subprocess_cmd('E:/propeller/qprop/bin/qprop E:/propeller/qprop/cases/input E:/propeller/qprop/cases/motor {} {} > {}'.format(self.V_in,self.rpm, name))
        output = subprocess_cmd('E:/propeller/qprop/bin/qprop E:/propeller/qprop/cases/input E:/propeller/qprop/cases/motor {} {}'.format(self.V_in,self.rpm))
        out = output.split()
        try:
            self.k_v = float(out[21])
            self.rho = float(out[28])
            self.miu = float(out[33])
            self.v = float(out[83])
            self.rpm = float(out[84])
            self.T = float(out[86])
            self.Q = float(out[87])
            self.P = float(out[88])
            self.Volts = float(out[89])
            self.Amps = float(out[90])
            results = {'v': self.v, 'rpm': self.rpm, 'T':self.T, 'Q':self.Q, 'P':self.P, 'Volts':self.Volts, 'Amps':self.Amps, 'g/W':self.T/9.81*1000/self.P}
            self.results = pd.DataFrame.from_dict(results,orient='index',columns=['value'])
        except IndexError:
            # this is case for bad design
            self.P = 1e5
            self.Q = 100
            self.T = 0
        except ValueError:
            # this is case for bad design 2
            self.P = 1e5
            self.Q = 100
            self.T = 0            
        
        
    def plots(self):
        plt.figure(dpi = 200)
        plt.plot(self.rs, self.chords, 'k-o')
        plt.xlabel('radius (m)')
        plt.ylabel('chord (m)')
        plt.grid(True)
#        plt.axis('equal')
        plt.xlim([-self.D/2*.1, self.D/2*1.1])
        plt.tight_layout()
        
        plt.figure(dpi = 200)
        plt.plot(self.rs, self.betas, 'k-o')
        plt.xlabel('radius (m)')
        plt.ylabel('beta (deg)')
        plt.grid(True)
#        plt.axis('equal')
        plt.xlim([-self.D/2*.1,self.D/2*1.1])
        plt.tight_layout()
        
    def plotTop(self):
        
        fig, ax = plt.subplots(1, 1, dpi = 150)
        
        circle = plt.Circle((0, 0), 0.05, color='k', fill=False)
        ax.add_artist(circle)
        ax.plot([0,self.rs[-1]], [0, 0], 'k--')
        ax.fill_between(self.rs, self.chords*np.cos(np.radians(self.betas)) *.75, color = '.5')
        ax.fill_between(self.rs, self.chords*np.cos(np.radians(self.betas)) *-.25, color = '.4')
        ax.axis('equal')
        ax.grid(True, 'both', linestyle = '--', linewidth = .2)
        ax.set_xlim(0, self.rs[-1]+.1)
        ax.set_ylim(np.max(self.chords)*-.25-.05, np.max(self.chords)*.75+.05)
        ax.set_xlabel('$r\,[m]$')
        ax.set_ylabel('$[m]$')
        ax.set_title('top view')
        fig.tight_layout()
        plt.show()
        
    def CreateMotorFile(self):
        """ created simplified motor data file  """
        option = 1
        if option ==1:
            with open(self.motorname, 'w') as file:
                file.write('\n')
                file.write('Motor name\n')
                file.write('\n')
                file.write(' 1     ! motor type (brushed DC)\n')
                file.write('\n')
                file.write(' 2     ! Rmotor (Ohms)\n')
                file.write(' 2     ! Io     (Amps)\n')
                file.write(' 24   ! Kv     (rpm/Volt)\n')
        else:
            with open(self.motorname, 'w') as file:
                file.write('test motor\n')
                file.write('\n')
                file.write(' 2        ! motor type (brushed DC, higher-order model)\n')
                file.write('\n')
                file.write(' 2     ! R0     (Ohms)      \n')     
                file.write(' 0.5    ! Io0    (Amps) \n')   
                file.write(' 210   ! Kv     (rpm/Volt)\n')
                file.write(' 280   ! Kq     (Amps/N-m)*30/pi\n')
                file.write(' 1.0E-5   ! tau    (s)\n')
                file.write(' 5.7E-5   ! Io1    (Amp-s)\n')
                file.write(' 4.0E-8   ! Io2    (Amp-s^2)\n')
                file.write(' 0.052    ! R2     (Ohms/Amps^2)\n')
        
        
    def CreateDataFile(self):
        """ generate input file with global airfoil chracteristic """
        with open(self.filename, 'w') as file:
            file.write(' \n')
            file.write('Test case\n')
            file.write('\n')
            file.write(' {}           ! Nblades\n'.format(self.N))
            file.write('\n')
            file.write(' 0.42  4.93   ! CL0     CL_a\n')
            file.write(' -0.1  .9   ! CLmin   CLmax\n')
            file.write('\n')
            file.write(' 0.0102  0.01295  0.0188  0.211  !  CD0    CD2u  CD2l    CLCD0\n')
            file.write(' 720000   .0       !  REref  REexp\n') # reynolds correction for large turbulent props goes to -.2, smaller props feel better with .7
            file.write('\n')
            file.write('\n')
            file.write(' 1.0  1.0   1.0  !  Rfac   Cfac   Bfac  \n')
            file.write(' 0.0     0.0      0.0  !  Radd   Cadd   Badd  \n')
            file.write('\n')
            file.write('#  r    chord    beta\n')
            for i in range(len(self.r)):                       
                file.write(' {}    {}    {}\n'.format(self.rs[i] , self.chords[i], self.betas[i] ))

    def CreateDataFile2(self):
        """ generate input file with locla airfoil characteristics based on Re """
        with open(self.filename, 'w') as file:
            file.write(' \n')
            file.write('Test case\n')
            file.write('\n')
            file.write(' 2           ! Nblades\n')
            file.write('\n')
            file.write(' 0.42  4.93   ! CL0     CL_a\n')  # these are global coefficients
            file.write(' -0.1  .7   ! CLmin   CLmax\n')   # ''
            file.write('\n')
            file.write(' 0.0102  0.01295  0.0188  0.211  !  CD0    CD2u  CD2l    CLCD0\n')
            file.write(' 720000   .0       !  REref  REexp\n') # reynolds correction for large turbulent props goes to -.2, smaller props feel better with .7
            file.write('\n')
            file.write('\n')
            file.write(' 1.0  1.0   1.0  !  Rfac   Cfac   Bfac  \n')
            file.write(' 0.0     0.0      0.0  !  Radd   Cadd   Badd  \n')
            file.write('\n')
            file.write('#  r    chord    beta\n')
            for i in range(self.N_stations):       
                #  r  chord  beta [  CL0  CL_a   CLmin CLmax  CD0   CD2u   CD2l   CLCD0  REref  REexp ]
                file.write(' {}    {}    {} {} {} {} {} {} {} {} {}\n'.format(self.rs[i] , self.chords[i], self.betas[i],
                                                   self.Cl0[i], self.Cla[i], self.Clmin[i], self.Clmax[i], 
                                                   self.Cd0[i], self.Cd2u[i], self.Cd2l[i], self.ClCd0[i], self.Reref[i]   ))



class Rotor2():
    """
    rotor from set of parameters defining chord and pitch distributions
    
    
    rpm = 1400
    x= [0.17745232142327558, 0.18625228204974623, 0.10804016368930469, 0.10197388288358913, 0.053621652353998905, 17.264483310376505, 16.55472426179523, 13.152857247471523, 6.964947458526998, 2.6881663537721785]
    pts = 12
    rot= Rotor2(x, 0.95, points = pts)
    betas = rot.betas
    chords = rot.chords
    rs = rot.Rs
    ts = np.linspace(0.024,.005,pts)
    den = 1
    rotor = Rotor(rpm = rpm, betas = betas[::den], rs = rs[::den], chords = chords[::den], ts = ts[::den])
    rotor.plotTop()

    """
    
    
    def __init__(self, x, r_max, typ = 1, points = 24):
        """ generate rotor from params from optimization (or in more general: with b-spline y-knots positions)
        
        if typ==1: require 10 parameters in x:
            - 5 for chord and
            - 5 for twist distributions
        with Rotor2 icem model can be generated
        
        """
        
        if typ == 1:
            num = 5
        else:
            num = 4
        
        # define knots for defining spline
        chord = x[:num]
        # const chord:
        
        beta = x[num:2*num]
        
        # got knots parameters, evaluate chord and beta distributions
        knots1 = np.zeros([num,2])
        knots1[:,0] = np.linspace(.2,r_max,num)
        
        knots2 = np.zeros([num,2])
        knots2[:,0] = np.linspace(.2,r_max,num)
        
        if typ == 1 or typ == 2:
            knots1[:,1] = chord
            knots2[:,1] = beta
        c1 = curve(knots1, points)
        c2 = curve(knots2, points)

        self.Rs = c1.curve[0]
        self.chords = c1.curve[1]
        self.betas = c2.curve[1]



class Rotor3():
    """ generate rotor with output file
    
        rpm = 1400
        rot = Rotor3(typ = 'type1-new',r = 0.95, rpm = rpm)
        betas = rot.betas
        chords = rot.chords
        rs = rot.Rs
        ts = np.linspace(0.024,.005,25)
        den = 1
        rotor = Rotor(rpm = rpm, betas = betas[::den], rs = rs[::den], chords = chords[::den], ts = ts[::den])
        rotor.plotTop()
            
            
            
    """
    def __init__(self, typ = 'type3', r = 0.85, rpm = 1498, nameprefix = 'type1-hybrid-'):
        
        filename =  r'E:/propeller/qprop/cases/results/'+ typ+ '/{}{}-{}.txt'.format(nameprefix, r, int(rpm))
        
        X = np.loadtxt(filename, skiprows = 20)

        self.Rs = X[:,0]
        self.chords = X[:,1]
        self.betas = X[:,2]











class curve:
    """
    cur = curve( np.array([[150,250],[150,300],[100, 200], [50, 150], [30, 100]]), 100)
    cur.savetxt()
    
    cur = curve( np.array([[200,250],[400,300],[600, 200], [800, 150], [1000, 100]]), 100)
    cur.savetxt()
    
    cur = curve( np.array([[200,150],[400,150],[600, 100], [800, 50], [1000, 30]]), 100)
    cur.savetxt()
    
    
    
    self.curve defines bspline from input knots
    
    """
    def __init__(self, knots, points=10):
        """ finds bspline curve for specified knots
        inputs: knots coords as numpy array  [  [x0,y0],[x1,y1],...,[xn-1,yn-1]  ]
        attribute: coords of points on spline, number defined in __init__  [ [x0,x1,x2,...,xn-1], [y0,y1,y2,...,yn-1] ]
        """
        self.points = points
        self.x=knots[:,0]
        self.y=knots[:,1]
                
        l=len(self.x)  
        
        t=np.linspace(0,1,l-2,endpoint=True)
        t=np.append([0,0,0],t)
        t=np.append(t,[1,1,1])
        tck=[t,[self.x,self.y], 3]
        u3=np.linspace(0,1,(max(l*2,points)),endpoint=True)
        self.curve = interpolate.splev(u3, tck, 0)
        
        
    def savetxt(self):
        xs =  self.curve[0].reshape(self.points, 1)
        ys =  self.curve[1].reshape(self.points, 1)
        
        np.savetxt('curvePoints.txt', np.append(xs, ys, axis = 1))
        
        xs = self.x.reshape(5, 1)
        ys = self.y.reshape(5, 1)
        np.savetxt('curveKnots.txt', np.append(xs, ys, axis = 1))
        
#if __name__ == '__main__':
#    rotor = Rotor(1.5)

    
#subprocess_cmd('pushd E:/propeller/qprop/runs&E:/propeller/qprop/bin/qprop cam6x3 s400-6v-dd 0 0 8')
    

    
    
def generate(airfoil, name):
    """ 
    generates airofil data for specified reynolds number range using XFOIL 
    
    input accepts airfoil object of class Airfoil
    and name of saved file
    
    """
    rel = np.arange(1e5, 1.4e6, 1e4)    
    reu = rel + 1e4
    n = len(rel)
    
    rel = rel.reshape((n, 1))
    reu = reu.reshape((n, 1))
    
    
    cl0 = np.zeros((n,1))
    cla = np.zeros((n,1))
    cd0 = np.zeros((n,1))
    clcd0 = np.zeros((n,1))
    cd2u = np.zeros((n,1))
    cd2l = np.zeros((n,1))

    
    for i in range(n):
        cl0[i], cla[i], cd0[i], clcd0[i], cd2u[i], cd2l[i] = airfoil.qpropData(m = 0, re = (rel[i,0]+reu[i,0])/2)
    
    print(cl0)
    
    X = np.concatenate( (rel, reu, cl0, cla, cd0, clcd0, cd2u, cd2l ) , axis = 1)
    
    print(X)
    np.savetxt(name, X)
    
def readQPolar( re, name = None, X = None):
    """ reads properties from matrix or text file
        aim to provide matrix, especially in optimization cases (efficiency)
    """
    if X is None:
        X = np.loadtxt(name)
    for rey in X:
        if re >= rey[0] and re < rey[1]:
            return rey[2:]
        
def generator():
    """
    generates airofil data for defined set of airfoils
    consider paths
    
    """
    for name in ['mh114','mh115','mh117'] :
        airfoil = Airfoil('XFOIL', r'E:/propeller/mh_airofils/{}/{}xfoil.txt'.format(name, name) )
        generate(airfoil, r'E:\propeller\qprop\airfoildata\{}.txt'.format(name)) 

    
    
