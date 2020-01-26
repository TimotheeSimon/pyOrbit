import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import os



from sys import path
path.append('D:\\Users\\s606482\\Desktop\\pyOrbit-master\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import planetary_data as pd
import tools as t


dt = 1
tspan = dt*3600*24*30
cb = pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2']=True
    op = OP(t.tle2coes('TLE\\AISSAT2.txt'), tspan, dt, coes=True, deg=False, perts=perts)
    op.plot_3d_animation()
