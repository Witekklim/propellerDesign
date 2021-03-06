# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:57:29 2019
@author: Witold Klimczyk

# ICEM
foil = Airfoil(filein = r'E:\propeller\mh_airofils\mh117/mh117.txt', t = 0.001, chord = 0.2)
foil.runFluent(15,.2,1)#

# XFOIL
foil2 = Airfoil(ftype = 'XFOIL', filein = r'E:\AIRFOIL\airfoils/naca0012.txt', t = 0.001, chord = 0.2)

# x,y


X = pd.read_csv(f'http://airfoiltools.com/airfoil/seligdatfile?airfoil={foilname}-il')
X.to_csv(r'E:\AIRFOIL\temp.csv', header = False, index = False)
X = np.loadtxt(r'E:\AIRFOIL\temp.csv')
foil = Airfoil( 'XFOIL',  r'E:\AIRFOIL\temp.csv')





"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
import subprocess
import os
from subprocess import run, PIPE
import gc
import pandas as pd
from urllib.error import HTTPError

from fluentScheme import generateScheme
from icemScheme import generateICEMScheme


class Airfoil():
    def __init__(self, ftype = 'ICEM', filein = None, x = None, y = None,  T_req = None, camber = None,
                 chord = None, beta = None, z = 0, fileoutICEM = None, t = 0, dx = 0, dy = 0, split = False, origin = 0,
                 camb = False, r_LE = None, verbose = False, workingdir = r'E:\AIRFOIL', xfoildir = None):
        """
        inputs:
            - ftype/name: 'ICEM' / 'XFOIL' / 'XY', specifies type of airofil input data or its name to download from airfoiltools
            - filein: '.txt' file with points coords (can be non-txt)
            - chord: dimensionless chord
            - beta: originally used for propeller pitch, for wing stands for twist
            - z: specifies third coordinae used for 3d wing stacking
            - fileoutICEM: full path and name for ICEM output file format, no extension (only name)
            - TEcut: specifies location of vertical cut
            - t: float: te thickness
            - T_req: maximum thickness to match required absolute thickness
            - origin: float: used to keep particular airfoil poitn in center, e.g. origin = .25 keeps quarter chord in center
            - camb: True/False: if we want to scael camber with thickness
            - workingdir: specify if other than current
            - xfoildir: contains xfoil.exe and uses this directory to save .txt files, if not given assumes it is in folder XFOIL under the same directory as current working folder
        attributes:
            - x: x-coords
            - y: y-coords
            - z: z-coords        
        """
        gc.collect()
        self.camber = camber
        self.chord = chord
        self.z = z
        self.filein = filein

        self.workingdir = workingdir if workingdir != None else os.getcwd().strip('\\python')
        print('workingdir {}'.format(self.workingdir))
        self.xfoildir = self.workingdir + '/XFOIL/'
        self.filebuffer = self.xfoildir + '/xfoilairfoil.txt'
        self.filecp = self.xfoildir + '/xfoilairfoilcp.txt'
        self.fileCptex = self.xfoildir + '/xfoilairfoilcptex.txt'
        self.camber_t = self.xfoildir + '/camber_t.txt'
        self.xfoilpath = self.xfoildir + '/xfoil.exe'

        # directories to check before analysis
        self.fileFig = self.workingdir + '/saved_plots/airfoil'
        self.meshin = self.workingdir + '/mesh/fluent.msh'
        self.meshDir = self.workingdir + '/mesh/'
        self.fileoutICEM = self.workingdir+'/mesh/foilICEM' if fileoutICEM is None else fileoutICEM
        self.fluentdir = self.workingdir + '/fluent/'
        self.ftype = ftype
        self.verbose = verbose
        
        self.camber = None
        self.thickness = None
        self.split = False
        self.t = t
        
        if ftype == 'ICEM':
            self.readICEM()
        elif ftype == 'XFOIL':
            self.readXFOIL()
            self.saveXFOIL()
        elif ftype == 'XY':
            self.x = x
            self.y = y
        else:
            try:
                X = pd.read_csv(f'http://airfoiltools.com/airfoil/seligdatfile?airfoil={ftype}-il')
                X.to_csv(r'E:\AIRFOIL\temp.csv', header = False, index = False)
                X = np.loadtxt(r'E:\AIRFOIL\temp.csv')
                self.x = X[:,0]
                self.y = X[:,1]
                self.z = self.z
            except HTTPError:
                print('error reading airofil from web')
                return None


        
        # chord scaling
        if chord is None:
            self.chord = np.max(self.x) - np.min(self.x)
            if self.verbose:
                print('evaluated chord is {:.2f}'.format(self.chord))
        else:
            self.chord = np.max(self.x) - np.min(self.x)
            self.scale_XFOIL(chord/np.max(self.x))
            self.saveXFOIL()
            if self.verbose:
                print('scaled airfoil to desired chord {:.3f}'.format(self.chord))

        self.x1 = None
        self.x2 = None
        self.y1 = None
        self.y2 = None
        self.z1 = z
        self.z2 = z
        
        # imposing required thickness
        if T_req is not None:
            self.thicken(T_req, camb)
        
        # cut TE
        if t > 0:
            self.cutTE_XFOIL(t, r = .5)
    
        if r_LE is not None:
            r_LE_current = self.LEradius()[0]
            if r_LE > r_LE_current:
                print('modifying LE radius')
                factor = r_LE / r_LE_current
                self.modify_XFOIL(1,1,factor)
                print('LE factor = {:.1f}'.format(factor))
        
        # twisting airfoil to required beta
        if beta is not None:
            self.rotate_XFOIL(beta, origin)
        
        # translating airofil to match required origin
        self.translate_XFOIL(dx - origin * self.chord, dy)
        
        # setting split after all airofil modifications
        self.split = split
        if self.split:
            self.splitCurve()
    
    def readICEM(self):
        X = np.loadtxt(self.filein, delimiter = '\t', skiprows = 1)
        self.x = X[:,0]
        self.y = X[:,1]
        self.z = self.z
    
    def readXFOIL(self, file = None):
        if file is None:
            X = np.loadtxt(self.filein, skiprows = 1)
        else:
            X = np.loadtxt(file, skiprows = 1)
        self.x = X[:,0]
        self.y = X[:,1]
        self.z = self.z
        self.filein = self.filebuffer
        
    def saveXFOIL(self):
        """  saves airfoil coords to .txt file with specified path 
        """
        # close trailing edge
        # save coords to file
        if not self.split:
            with open(self.filebuffer, "w") as text_file:
                print("airfoil", file = text_file)
                for i in range(len(self.x)):
                    print("  {} {}".format(self.x[i], self.y[i]), file=text_file)
        else:
            with open(self.filebuffer, "w") as text_file:
                print("airfoil", file = text_file)
                for i in range(len(self.x2)-1,0,-1):
                    print("  {} {}".format(self.x2[i], self.y2[i]), file=text_file)
                for i in range(len(self.x1)):
                    print("  {} {}".format(self.x1[i], self.y1[i]), file=text_file)

###
###   ===================        GEOMETRY SECTION           ==========================
###
                    
    def cutTE_XFOIL(self, t = .005, r = 0.5):
        """ modifies airfoil using xfoil to maintain camber
        
            t: thickness
            r: blending radius
        """
        self.saveXFOIL()
        airfoilIN = self.filebuffer 
        airfoilOUT = self.filebuffer
        command = 'load ' + airfoilIN + '\npane\ngdes\ntgap '+ '{} {}'.format(t,r) + '\n\npane\n\nsave '+airfoilOUT+'\ny\n\nquit\n'
        run([self.xfoilpath], stdout=PIPE,  input=command, encoding='ascii', shell = False)
        if self.verbose:
            print('succesfully modified TE using xfoil')
        self.readXFOIL(airfoilOUT)
        
    def modify_XFOIL(self, thicken = 1, camber = 1, LE_radius=1):
        """ modifies airfoil using xfoil to scale"
        thickness and camber distribution
            values below 1 decrease, above 1 increase scale (1 is no change)
        """
        self.saveXFOIL()
        airfoilIN = self.filebuffer 
        airfoilOUT = self.filebuffer 
        command = 'load ' + airfoilIN + '\npane\ngdes\ntfac '+ '{} {}'.format(thicken, camber) +'\nlera {} {}'.format(LE_radius, .2)+  '\n\npane\n\nsave '+airfoilOUT+'\ny\n\nquit\n'

        p = run([self.xfoilpath], stdout=PIPE,  input=command, encoding='ascii', shell = False)
        if self.verbose:
            print('modified thickness scaled by {}'.format(thicken))
        self.readXFOIL(airfoilOUT)    
        
    def thicken(self, req_T, camb = False):
        """ modifies thickness to required value of maximum thickness
            can also modify camber of airfoil
        """
        self.findCamberThickness()
        factor = req_T/(self.t_max*self.chord)
        print(f'{factor}')
        if camb==True:
            camb = factor
            self.modify_XFOIL(thicken = factor, camber = camb)
        else:
            camb = 1
            self.modify_XFOIL(thicken = factor, camber = camb)
        if self.verbose:
            print('modified thickness to desired value, i.e. {:.3f}, by a factor of {:.2f}'.format(req_T, factor))
            
        
    def scale_XFOIL(self, factor = 1):
        """ scales airfoil using xfoil 
        """
        print('chord before modification: {:.3f}'.format(self.chord))
        self.saveXFOIL()
        airfoilIN = self.filebuffer 
        airfoilOUT = self.filebuffer 
        command = 'load ' + airfoilIN + '\npane\ngdes\nscal '+ '{}'.format(factor) + '\n\npane\n\nsave '+airfoilOUT+'\ny\n\nquit\n'

        p = run([self.xfoilpath], stdout=PIPE,  input=command, encoding='ascii', shell = False)
        
        if self.verbose:
            print('modified chord by factor {}'.format(factor))
        self.readXFOIL(airfoilOUT)    
        self.chord *= factor
        print('chord after modification: {:.3f}'.format(self.chord))
        
    def translate_XFOIL(self, dx = 0, dy = 0):
        """ translates airfoil by specified dx and dy
        """
        self.saveXFOIL()
        airfoilIN = self.filebuffer 
        airfoilOUT = self.filebuffer 
        command = 'load ' + airfoilIN + '\npane\ngdes\ntran '+ '{} {}'.format(dx, dy) + '\n\npane\n\nsave '+airfoilOUT+'\ny\n\nquit\n'

        p = run([self.xfoilpath], stdout=PIPE,  input=command, encoding='ascii', shell = False)
        if self.verbose:
            print('airfoil translated by {:.3f} in x and {:.3f} in y'.format(dx, dy))
        self.readXFOIL(airfoilOUT)   
        
    def rotate_XFOIL(self, angle = 0, origin  = 0):
        """ rotates airfoil using xfoil by specified angle in degrees, around (0,0), positive angle moves TE down
        """
        if origin is not 0:
            self.translate_XFOIL(dx = -origin*self.chord)
        
        
        self.saveXFOIL()
        airfoilIN = self.filebuffer 
        airfoilOUT = self.filebuffer 
        command = 'load ' + airfoilIN + '\npane\ngdes\nadeg '+ '{}'.format(angle) + '\n\npane\n\nsave '+airfoilOUT+'\ny\n\nquit\n'

        p = run([self.xfoilpath], stdout=PIPE,  input=command, encoding='ascii', shell = False)
        if self.verbose:
            print('airfoil rotated by {:.2f}'.format(angle))
        self.readXFOIL(airfoilOUT)
        
        if origin is not 0:
            self.translate_XFOIL(dx = origin*self.chord )       
            
    def findCamberThickness(self, plot = False, tex = False, name = ''):
        """ finds camber and thickness distributions usign xfoil  """

        self.saveXFOIL()
        airfoilIN = self.filebuffer 
        airfoilOUT = self.filebuffer 
        command = 'load ' + airfoilIN + '\npane\ngdes\ntcpl\ncamb\nwrtc\n{}\n\n\nquit\n'.format(self.camber_t)

        p = run([self.xfoilpath], stdout=PIPE,  input=command, encoding='ascii', shell = False)
        if p.returncode ==2:
            if self.verbose:
                print('found camber and thickness distributions')
            
        X = np.loadtxt(self.camber_t, skiprows = 1)
        self.camber = X[:,:2]
        self.thickness = X[:,2:]
        self.readXFOIL(airfoilOUT)
        self.t_max = 2* np.max(self.thickness[:,1])
        
        
        if plot:        
            plt.figure(figsize = (6,2),dpi = 200)
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')

            plt.plot(self.camber[:,0], self.camber[:,1], 'k-',linewidth = 1.2, label = 'camber')
            plt.plot(self.thickness[:,0], self.thickness[:,1], 'k--',linewidth = 1.2, label = 'thickness')
            plt.plot(self.thickness[:,0], -self.thickness[:,1], 'k--',linewidth = 1.2)
            plt.xlabel(r'$x/c$')
            plt.ylabel(r'$y/c$',fontsize=12)
            plt.title(r"{}".format('camber and thickness distributions'), fontsize=12)
            #plt.subplots_adjust(top=0.8)
            plt.axis('equal')
            plt.legend()
            plt.tight_layout()
            plt.grid('major', linewidth = .2)
            plt.savefig(self.fileFig+'ct', dpi = 1000)
            plt.show()
            
        if tex:
            camberdir = self.workingdir + r'\wing3d\tex-plots\{}camber.txt'.format(name)
            thicknessdir = self.workingdir + r'\wing3d\tex-plots\{}thickness.txt'.format(name)
            np.savetxt(camberdir, self.camber)
            np.savetxt(thicknessdir, self.thickness)

    
    def t_x(self, x=None):
        """ finds thickness at specified x-position, x is x/c (i.e. between 0-1)
            
        self.t_x(0.5) returns thickness at x/c = 0.5
        if no argument passed, returns max thickness
        """
        self.findCamberThickness()
        i = 0
        if x is None:
            return self.t_max
        
        for i in range(len(self.thickness[:,0])):
            if self.thickness[i,0] > x :
                return 2*self.thickness[i,1]
        if self.verbose:
            print('invalid argument')
        
        
    
            
    def LEradius(self, plot = False, dpi = 500, saveFig = False):     
        """ method to find leading edge radius 
        
            buids many circles, each from 3 points from leading edge region
            lowest radius circle is chosen as le radius

            allows to plot le region to investigate le radius
        """
        def findCircle(P1, P2, P3):
            import sympy as sym
            
            a, b, r2 = sym.symbols('a, b, r2')
            e1 = sym.Eq((P1[0]-a)**2+(P1[1]-b)**2, r2**2)
            e2 = sym.Eq((P2[0]-a)**2+(P2[1]-b)**2, r2**2)
            e3 = sym.Eq((P3[0]-a)**2+(P3[1]-b)**2, r2**2)
            
            solution = sym.solve([e1, e2, e3], (a, b, r2))
            r = float(np.abs(solution[0][2]))
            x = float(np.abs(solution[0][0]))
            y = float(np.abs(solution[0][1]))
            return x,y,r
        
        
        i = np.where(self.x == min(self.x))[0][0]
        # find several circles around LE
        r = 1
        j = 1
        k = 1
        while j<5:
            while k<5:
                x_temp,y_temp,r_temp = findCircle( [self.x[i-j], self.y[i-j]], [self.x[i], self.y[i]], [self.x[i+k], self.y[i+k]]  )
                if r_temp<r:
                    r = r_temp
                    x = x_temp
                    y = y_temp
                k+=1
            j+=1

        if plot:
            an = np.linspace(0, 2*np.pi, 100)
            
            plt.figure(dpi = dpi)
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.plot(self.x, self.y, 'ko-', linewidth = 1.4)
            plt.plot([x],[y],'ro')
            plt.plot(r*np.cos(an)+x, r*np.sin(an)+y, 'r-', linewidth = 1.4)
            plt.title(r"{}".format('leading edge radius close up'), fontsize=12)
            
            plt.axis('equal')
            plt.ylim(-r, r*3.5)

            if saveFig:
                plt.savefig(self.fileFig, dpi = 1000)
            plt.show()
                
            fig, ax = plt.subplots(dpi = 500)
            ax.plot(self.x, self.y, 'ko-', linewidth = 1.4)
            ax.plot([x],[y],'ro')
            ax.plot(r*np.cos(an)+x, r*np.sin(an)+y, 'r-', linewidth = 1.4)
            ax.set_xlim(-r, r*3.5)
            ax.set_ylim(-r*2, r*2)
            ax.set_title('mh117: R=2')
            ax.set_aspect(1.0)
            ax.grid(which='major', linewidth = 0.2)
            plt.show()
            
        return r, x , y
        

    
    def saveICEM(self, airfoilfile = None):
        """ saves points in icem format, either as a single curve of splits to upper and lower (recommended) """
        
        if self.y[1]>self.y[-1]:
            self.x = np.flip(self.x, axis = 0)
            self.y = np.flip(self.y, axis = 0)
        
        if airfoilfile is not None:
            self.fileoutICEM = airfoilfile
            
        if not self.split:
            self.zs = np.ones(len(self.x))*self.z
            self.fileoutICEM += '.txt'
            with open( self.fileoutICEM, 'w') as f:
                f.write('{}\t{}\n'.format(len(self.x), 1))
                for i in range(len(self.x)):
                    f.write('{}\t{}\t{}\n'.format(self.x[i]*1000, self.y[i]*1000, self.zs[i]*1000) )
        else:
            self.z1 = np.ones(len(self.x1))*self.z
            self.z2 = np.ones(len(self.x2))*self.z
            
            with open( self.fileoutICEM + '.0.txt', 'w') as f:
                f.write('{}\t{}\n'.format(len(self.x1), 1))
                for i in range(len(self.x1)):
                    f.write('{}\t{}\t{}\n'.format(self.x1[i]*1000, self.z1[i]*1000, self.y1[i]*1000) )
                    
            with open( self.fileoutICEM + '.1.txt', 'w') as f:
                f.write('{}\t{}\n'.format(len(self.x2), 1))
                for i in range(len(self.x2)):
                    f.write('{}\t{}\t{}\n'.format(self.x2[i]*1000, self.z2[i]*1000, self.y2[i]*1000) ) 
                    
    def saveSW(self, airfoilfile):
        """ saves points in sw format, either as a single curve of splits to upper and lower (recommended) """
            
        if not self.split:
            self.zs = np.ones(len(self.x))*self.z
            airfoilfile += '.txt'
            with open( airfoilfile, 'w') as f:
                for i in range(len(self.x)):
                    f.write('{}\t{}\t{}\n'.format(self.x[i]*1000, self.zs[i]*1000, self.y[i]*1000) )
        else:
            self.z1 = np.ones(len(self.x1))*self.z
            self.z2 = np.ones(len(self.x2))*self.z
            
            with open( airfoilfile + '.0.txt', 'w') as f:
                for i in range(len(self.x1)):
                    f.write('{}\t{}\t{}\n'.format(self.x1[i]*1000, self.z1[i]*1000, self.y1[i]*1000) )
                    
            with open( airfoilfile + '.1.txt', 'w') as f:
                for i in range(len(self.x2)):
                    f.write('{}\t{}\t{}\n'.format(self.x2[i]*1000, self.z2[i]*1000, self.y2[i]*1000) ) 
                    
                       
###
###   ===================        ANALYSIS SECTION           ==========================
###
   
        
    def runXFOIL(self, cl=.2, alfa = None, re=1e6, m =.2, n_crit = 6, iters = 500, cp = False):
        self.saveXFOIL()
        airfoilIN = self.filebuffer 

        if alfa is None:
            S = cl
            s = 'cl'
            if self.verbose:
                print('running XFOIL for: cl={}'.format(cl))
        else:
            S = alfa
            s = 'a'
            if self.verbose:
                print('running XFOIL for: aoa={}'.format(alfa))
        if not cp:
            commands = 'load ' + airfoilIN + '\npane\noper\nvpar\nn {}\n\nvisc {}'.format(n_crit, re) + '\niter '+str(iters)+'\n{} {}'.format(s, S) + '\n\nquit\n'
            p = run([self.xfoilpath], stdout=PIPE,
                    input=commands, encoding='ascii')
        
        else:
            commands = 'load ' + airfoilIN + '\npane\noper\nvpar\nn {}\n\nvisc {}'.format(n_crit, re) + '\niter '+str(iters)+'\n{} {} '.format(s, S) + '\ncpwr\n{}\n\nquit\n'.format(self.filecp)
            p = run([self.xfoilpath], stdout=PIPE,
                    input=commands, encoding='ascii')
            return 0

        try:
            alfa = float(p.stdout[-130:-118])
            Cl = float(p.stdout[-112:-106])
            Cd = float(p.stdout[-78:-69])
            Cm = float(p.stdout[-94:-86])
            
            print(alfa,Cl,Cd,Cm)
        except ValueError:
            if self.verbose:
                print('error running xfoil, try slighlty different cl/alpha') # the reason is xfoil may not converge for this particular condition but in general it converges
            if alfa is None:
                alfa, Cd, Cm, Cl = self.runXFOIL(cl = 1.01*cl, re = re, m = m, n_crit = n_crit, iters = iters)
            else:
                alfa, Cd, Cm, Cl = self.runXFOIL(alfa = .01+alfa, re = re, m = m, n_crit = n_crit, iters = iters)
                
            #return 1, 1, 1, 1
            
        return alfa, Cd, Cm, Cl
    
    
    def runPolar(self, a0=-4, a1=8, re=1e6, m=.2, n_crit = 6, plot = False):
        alfas = np.zeros(a1-a0)
        cds = np.zeros(a1-a0)
        cls = np.zeros(a1-a0)
        cms = np.zeros(a1-a0)
        i=0
        for aoa in np.arange(a0,a1,1):
            alfas[i], cds[i], cms[i], cls[i] = self.runXFOIL(alfa = aoa, re= re, m = m, n_crit = n_crit)
            i+=1
            
        alfas = np.delete(alfas, np.where(cds == 1))
        print (alfas)
        
        if plot:
            plt.figure()
            plt.plot(alfas, cls, 'o-')
            plt.xlabel(r'$ \alpha [^\circ]$')
            plt.ylabel(r'$C_L$')
            plt.show()
            plt.figure()
            plt.plot(alfas, cds, 'o-')
            plt.xlabel(r'$ \alpha[^\circ]$')
            plt.ylabel(r'$C_D$')
            plt.show()
        
        return alfas, cds, cms, cls
    
    def plotCp(self, outputtex = False, dpi = 200, name = None, saveFig = False, airfoil= True, alfa = None):
        X = np.loadtxt(self.filecp, skiprows = 3)
        x = X[:,0]
        cp = X[:,2]
        if outputtex:
            np.savetxt(self.fileCptex, X)
        
        if name is None:
            if alfa is not None:
                name = '$C_p$ distribution at $\alpha = {}$'.format(alfa)
            else:
                name = '$C_p$ distribution'
        
        plt.figure(figsize = (6,4),dpi = dpi)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        plt.plot(x, -cp, 'k-',linewidth = 1)
        
        if airfoil:
            plt.plot(self.x/self.chord, self.y/self.chord*3-np.max(cp), 'k-',linewidth = 1)
        
        plt.xlabel(r'$x/c$',fontsize=12)
        plt.ylabel(r'$-C_p$',fontsize=12)
        plt.title(r"{}".format(name), fontsize=12)
        plt.subplots_adjust(top=0.8)
#        plt.axis('equal')
        plt.grid(which='major', linewidth = 0.2)
        plt.tight_layout()
#        plt.grid(True)
        if saveFig:
            plt.savefig(self.fileFig, dpi = 1000)
        plt.show()
        

    def runFluent(self, alfa, mach, chord, 
                  rho = 1.225, T = 300, viscosity = 1.78e-5,
                  name = 'airfoil', path = None, ID = 0,
                  mesh = 'o', y1 = 0.01, n_r = 120, n_le = 30, n_top = 120,
                  model = 'kw-sst', intermittency = False, lowre = False, polar = False,
                  onlymesh = False, onlyfluent = False, mshin = None, meshunits = 'mm',
                  tt = 1, farfieldnames = ['farfield'], outletnames = [], interiornames = ['int_fluid']
                  ):
        """
        chord used to scale mesh in fluent and use for coefficients
        if using auto o-mesh, generate airfoil with unit chord and scale mesh to required value
        
        
        
        static method: can be applied for given mesh, without airfoil initialization
        
        
        """
        if path is None:
            path = self.workingdir + r'\fluent'
        
        import time
        start = time.time()
        # begin with structured mesh generation
        # 
        import subprocess
        def subprocess_cmd(command):
            process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
            proc_stdout = process.communicate()[0].strip()
            #    print(proc_stdout)
            return proc_stdout
        
        if not onlyfluent:
            self.saveICEM(self.fileoutICEM)
            ICEMrun ='"C:\\Program Files\\ANSYS Inc\\v194\\icemcfd\\win64_amd\\bin\\icemcfd" -script'
            
                # pick mesh replay file to generate mesh
            if mesh == 'o':
                meshrpl = self.meshDir + 'omesh.rpl'
    #            ICEMscr = r'"E:\propeller\python\wing3d\rpl42.rpl"'
    #            ICEMscr = r'"E:\propeller\python\wing3d\omesh\omesh.rpl"'
                ICEMscr = f'"{meshrpl}"'
            elif mesh == 'unstructured':
                ICEMscr = r'"C:\Users\wk5521\Documents\ICEM\airfoil replays\mesh_output.rpl"'
            
            
            generateICEMScheme( y1 = y1, n_r = n_r, n_le = n_le, n_top = n_top, file = meshrpl)
            
            ICEM = ICEMrun + ' ' + ICEMscr
            subprocess_cmd(ICEM)
        
        # now having the mesh, run shceme generation, hence fluent
        if onlymesh:
            print('finished mesh')
            return 0
        
        
        fluentjournal = self.workingdir + '/fluent/journal.txt'
        
        casename = f'foil,{model},{alfa},{mach},{chord},{self.t}'
        if polar:
            casename = f'foil,{model},{mach},{chord},{self.t}'
        
        if lowre:
            casename+=',lowre'
        if intermittency:
            casename+= 'inter'
        
        meshin = mshin if mshin is not None else self.meshin
        
        generateScheme(filename = fluentjournal,
           casename = casename, 
           chord = chord, 
           viscosity = viscosity, 
           T=T,
           alfa = alfa,
           mach = mach,
           meshin = meshin,
           meshunits = meshunits,
           farfieldnames = farfieldnames,
           outletnames = outletnames,
           interiornames = interiornames,
           path = self.fluentdir,
           model = model,
           intermittency = intermittency,
           lowre = lowre,
           polar = polar, 
           tt =tt
           )
    
        FLUENTrun = '"C:\\Program Files\\ANSYS Inc\\v194\\fluent\\ntbin\\win64\\fluent.exe" 2d -t8 -wait -i'
        FLUENT  = FLUENTrun + ' '+ '"{}"'.format(fluentjournal)
        
        subprocess_cmd(FLUENT)
        end = time.time()
        
        showresult = False
        if showresult:
            result = np.loadtxt('{}/reports/{}.out'.format(self.fluentdir, casename), skiprows = 100)
            result = result[-10:]
            result = np.mean(result, axis = 0)
            lift = result[1]
            drag = result[2]
            moment = result[3]
            duration = end - start
            print('mesh size: {}, lift: {:.4f}, drag: {:.6f}, duration: {}'.format(2*(n_le+n_top)*n_r , lift , drag , duration))
            return 2*(n_le+n_top)*n_r , lift , drag

    def splitCurve(self):
        """ splits curve into two curves at leading edge by front-most point """
        i_min = np.where(self.x == np.amin(self.x))[0][0]
        self.split = True
        self.x1 = self.x[:i_min+1]
        self.y1 = self.y[:i_min+1]
        self.x2 = self.x[i_min:]
        self.y2 = self.y[i_min:]
        self.z1 = np.ones(len(self.x1))*self.z
        self.z2 = np.ones(len(self.x2))*self.z
        
    def qpropData(self, m, re, n = 12, n_crit = 5):
        """ this method finds coefficients required to define qprop input file            
        returns (cl0, clalfa, cd0, clcd0, cd2u, cd2l)
        """
#        collect some data for range of angles of attack
        alfas = np.zeros(n)
        cds = np.zeros(n)
        cms = np.zeros(n)
        cls = np.zeros(n)
        
        j = 0
        for i in range(-6, 6, 1):
            self.cutTE_XFOIL(t = 0.005, r = .3)
            alfas[j], cds[j], cms[j], cls[j] = self.runXFOIL(alfa = i, re = re, m = m, iters = 1000, n_crit = n_crit)
            j+=1
        
        cl0 = cls[6]
        clalfa = (cls[-1] - cls[4]) / np.radians(7)
        
        # now begin drag section
        from scipy.optimize import minimize, fmin
        cd0 = cds.min()

        for index in range(len(cds)):
            if cds[index] == cd0:
                break
        clcd0 = cls[index]
        
        def merit(x):
            # args = (cd0, )
            merit = np.abs(  cds[index + 1] + cds[index + 2] - (cd0 + x * (cls[index + 1] - clcd0 )**2 + cd0 + x * ( cls[index + 2] - clcd0 )**2 ) )
            return  merit
        
        result = fmin(merit,  .1)
        cd2u = result[0]
        
        def merit2(x, *args):
            return  np.abs(  cds[args[0] - 1] + cds[args[0] - 2] - (cd0 + x * (cls[args[0] - 1] - clcd0 )**2 + cd0 + x * ( cls[args[0] - 2] - clcd0 )**2 ) )
        
        result2 = minimize(merit2, .05, args = (index))
        cd2l = result2.x[0]

        print('cl0, clalfa, cd0, clcd0 = {:.3f} {:.3f} {:.3f} {:.3f}'.format( cl0, clalfa, cd0, clcd0))
        return cl0, clalfa, cd0, clcd0, cd2u, cd2l
        
    
    
    def plotAirfoil(self, name=None, saveFig = False, dpi = 200, tex = False , nametex = ''):
        
        if name is None:
            name = 'airfoil'
        
        plt.figure(figsize = (6,2),dpi = dpi)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        if self.split:
            plt.plot(self.x1/self.chord, self.y1/self.chord, 'k-',linewidth = 1.2)
            plt.plot(self.x2/self.chord, self.y2/self.chord, 'k-',linewidth = 1.2)
        else:
            plt.plot(self.x/self.chord, self.y/self.chord, 'k-',linewidth = 1.2)
        plt.xlabel(r'$x/c$',fontsize=12)
        plt.ylabel(r'$y/c$',fontsize=12)
        plt.title(r"{}".format(name), fontsize=12)
        plt.subplots_adjust(top=0.8)
        plt.axis('equal')
        plt.grid(which='major', linewidth = 0.2)
        plt.tight_layout()
        if saveFig:
            plt.savefig(self.fileFig, dpi = 1000)
        plt.show()
        
        if tex:
            X = np.append((self.x/self.chord).reshape(-1,1), (self.y/self.chord).reshape(-1,1), axis = 1 )
            savedir = self.workingdir + r'\wing3d\tex-plots\{}airfoil.txt'.format(nametex)
            np.savetxt(savedir, X)
        
