import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D



from sys import path
path.append('.\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import planetary_data as pd
import tools as t


tspan = 24*3600*10
dt = 60
cb = pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2']=True
    mass0 = 1000


    coes = np.array([cb['radius']+600, 0.01, 90.01, 0.0, 0.0, 0.0,[0,0,0]])
    op = OP(coes, tspan, dt, mass0, coes=True, deg=True, perts=perts)

    op.calculate_coes()


    op.plot_coes(hours=True, show_plot=True)
