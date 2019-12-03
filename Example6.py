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


tspan = 4*24*3600
dt = 100
cb = pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2']=True
    rv = np.array([-3670, -3870, 4400, 4.7, -7.4, 1])
    op = OP(rv, tspan, dt, coes=False, deg=False, perts=perts)

    # coes = np.array([cb['radius']+600, 0.1, 90.01, 0.0, 0.0, 0.0])
    #op = OP(coes, tspan, dt, coes=True, deg=True, perts=perts)

    op.calculate_coes()

    print((op.coes[-1,5]-op.coes[0,5])/tspan)
    print((op.coes[-1,4]-op.coes[0,4])/tspan)
    print(op.coes[-1,5]) 
    print(op.coes[-1,4]) 

    op.plot_coes(hours=True, show_plot=True)
    t.plot_3d(op.ys)
