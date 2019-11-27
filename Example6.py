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


tspan = 3600*24
dt = 10
cb = pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2']=True
    op = OP(t.tle2coes('Orbit Visualisation\\AISSAT2.txt'), tspan, dt, coes=True, deg=False, perts=perts)
    op.calculate_coes()
    op.plot_coes(hours=True, show_plot=True)

	