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



dt = 1
tspan = dt*60*24*36
cb = pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2']=False
    perts['thrust']=10.089
    perts['isp']=1650
    
    mass0 = 367

    coes = np.array([cb['radius']+1000, 0.1, 10, 0.0, 0.0, 0.0, [0,0,0]])
    op = OP(coes, tspan, dt, mass0, coes=True, deg=True, perts=perts)

    op.calculate_coes()
    op.plot_coes(hours=True)
    #op.plot_3d()
    op.plot_3d_animation()
    #op.plot_masses(hours=True)
