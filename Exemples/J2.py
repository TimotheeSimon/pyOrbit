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


tspan = 24*3600*2.
dt = 100.
cb = pd.earth

if __name__ == '__main__':
    perts = perturbations()
    perts['J2']=True

    #book example
    r0 = np.array([-2384.46,5729.01,3050.46])
    v0 = np.array([-7.36138,-2.98997,1.64354])
    state0=np.array(t.rv2coes(r0,v0,print_results=True, deg=True)+[[0,0,0]])

    #ISS example
    state0=t.tle2coes('TLE\\ISS.txt', deg=True)


    #no raan shift example
    state0=np.array([cb['radius']+600, 0.01, 90.01, 0.0, 0.0, 0.0,[0,0,0]])

    #no aop shift example
    state0=np.array([cb['radius']+600, 0.01, 63.435, 0.0, 0.0, 0.0,[0,0,0]])

    op = OP(state0, tspan, dt, coes=True, deg=True, perts=perts)
    op.plot_3d(show_plot=True)
    op.calculate_coes()
    op.plot_coes(hours=True, show_plot=True)