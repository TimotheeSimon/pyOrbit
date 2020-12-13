import numpy as np
from sys import path
path.append('.\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import perturbations
import planetary_data as pd
import tools as t


dt = 60
tspan = 3600*24
cb = pd.earth

if __name__ == '__main__':
    perts = perturbations()

    perts['raan_maneuvre']=True
    perts['thrust']=0.1
    perts['isp']=1650
    
    mass0 = 367

    coes = np.array([cb['radius']+1000.0, 0.01, 30.0, 10.0, 110.0, 300.0, [0,0,0]])
    op = OP(coes, tspan, dt, mass0, coes=True, deg=True, perts=perts)

    op.calculate_coes()
    op.plot_coes(hours=True, show_plot=True)
    #op.plot_3d_animation()
    #op.plot_masses(hours=True)
