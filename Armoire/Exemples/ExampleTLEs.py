import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D



from sys import path
path.append('D:\\Users\\s606482\\Desktop\\pyOrbit-dev\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
import planetary_data as pd
import tools as t

dt = 10
tspan = 24*600*dt

if __name__ == '__main__':
    cb = pd.earth
    op0 = OP(t.tle2coes('TLE\\luch.txt'), tspan, dt, coes=True, deg=False)
    op2 = OP(t.tle2coes('TLE\\AISSAT2.txt'), tspan, dt, coes=True, deg=False)
    #t.plot_n_orbits([op0.rs, op2.rs], labels=['Luch', 'AISSAT 2'], show_plot=True)
    t.plot_3d(op2.rs, show_plot=True)
