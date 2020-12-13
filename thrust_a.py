import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D



from sys import path
path.append('.\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import perturbations
import planetary_data as pd
import tools as t



dt = 2
tspan = 3600*24
cb = pd.earth

if __name__ == '__main__':
    perts = perturbations()
    perts['J2']=True
    perts['thrust']= 0.5
    perts['isp']=1650
    
    mass0 = 100
    
    # a,e,i,ta,aop,raan, date
    coes = np.array([cb['radius']+1000, 0.01, 0.1, 0.0, 0.0, 0.0, [0,0,0]])
    op = OP(coes, tspan, dt, mass0, coes=True, deg=True, perts=perts)

    op.calculate_coes()
    #op.plot_coes(hours=True, show_plot=True)
    op.plot_3d()
    #op.plot_3d_animation()
    #op.plot_masses(hours=True)
