import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import os



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
    perts['J2']=True
    coes = t.tle2coes('TLE\\AISSAT2.txt')
    op = OP(coes, tspan, dt, coes=True, deg=False, perts=perts)
    op.plot_3d()
