import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D



from sys import path
path.append('C:\\Users\\timsi\\Documents\\Visual Studio Code\\Orbit Visualisation\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import planetary_data as pd
import tools as t



dt = 600
tspan = 3600*24
cb = pd.earth

if __name__ == '__main__':

    N = 1

    perts = null_perts()
    perts['J2']=True
    perts['thrust']=0.09
    perts['isp']=1660
    perts['direction']=1
    
    t0 = 0
    mass0 = 1430
     
    
    M1_init = np.array([cb['radius']+700, 0.01, 90, -106, 170, 1.6, [0,0,0]])
    op = OP(M1_init, tspan, dt, mass0, t0, coes=True, deg=True, perts=perts)
    op.calculate_coes()
    op.plot_coes(days=True, show_plot=True)


