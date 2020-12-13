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


tspan = 24*3600
dt = 60
cb = pd.earth

if __name__ == '__main__':
    perts = perturbations()
    perts['J2']=True

    coes = np.array([cb['radius']+1000, 0.01, 0.01, 0.0, 0.0, 0.0,[0,0,0]])
    op = OP(coes, tspan, dt, coes=True, deg=True, perts=perts)

    op.calculate_coes()
    op.plot_coes(hours=True, show_plot=True)
