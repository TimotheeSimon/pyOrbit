import numpy as np
from sys import path
path.append('.\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import perturbations
import planetary_data as pd
import tools as t


# Pour cette manoeuvre, une poussée positive réduit l'eccentricité -> à verifier

dt = 60
tspan = 24*3600
cb = pd.earth

if __name__ == '__main__':
    
    perts = perturbations()
    perts['a_maneuvre']=False
    perts['ecc_maneuvre']=True
    perts['thrust']=0.089
    perts['isp']=1650
    
    mass0 = 367

    coes = np.array([cb['radius']+1800, 0.01, 90, 10, 0.0, 0.0, [0,0,0]])  #a, e_norm, i, ta, aop, raan
    op = OP(coes, tspan, dt, mass0, coes=True, deg=True, perts=perts)
    
    op.calculate_coes()
    op.plot_coes(hours=True, show_plot=True)
    op.plot_3d()
    op.plot_masses(hours=True, show_plot=True)

