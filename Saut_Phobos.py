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

dt = 60
tspan = 3600*24

if __name__ == '__main__':
    perts = perturbations()
    perts['thrust']= 0.5
    perts['isp']=1650
    perts['J2']=False
    
    mass0 = 65
    r_mag = pd.phobos['radius']
    
    state0 = np.array([r_mag, 0, 0, 0.1, 0, 0])
    op = OP(state0, tspan, dt, mass0, coes=False, cb=pd.phobos, perts=perts)
    op.calculate_coes()
    #op.plot_coes(hours=True, show_plot=True)
    #op.plot_3d()
    op.plot_3d_animation()
    #op.plot_masses(hours=True)
