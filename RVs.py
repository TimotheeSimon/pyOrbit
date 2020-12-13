import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D


from sys import path
path.append('.\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
import planetary_data as pd
import tools as t

dt = 60
tspan = 24*60*dt

if __name__ == '__main__':
    
    r_mag = pd.earth['radius']+35000
    v_mag= sqrt(pd.earth['mu']/r_mag)
    state0 = np.array([r_mag, 0, 0, 0, -v_mag*0.7, v_mag*0.2])

    r_mag = pd.earth['radius']+8000
    v_mag= sqrt(pd.earth['mu']/r_mag)*1.1
    state00 = np.array([r_mag, r_mag/2, 0, v_mag*0.4, -v_mag*0.8, -v_mag*0.2])

    op0 = OP(state0, tspan, dt)
    op00 = OP(state00, tspan, dt)

    op0.propagate_orbit()
    op00.propagate_orbit()

    t.plot_n_orbits([op0.rs, op00.rs], labels=['orbite A','orbite B'], show_plot=True)
    
