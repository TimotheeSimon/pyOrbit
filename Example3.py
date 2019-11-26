import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D


from sys import path
path.append('C:\\Users\\timsi\\Documents\\Visual Studio Code\\Orbit Visualisation\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
import planetary_data as pd
import tools as t

dt = 60
tspan = 24*60*dt

if __name__ == '__main__':
	cb = pd.earth
	# a,e,i,ta,aop,raan
	c0 = [cb['radius']+1202.0, 0.0001604, 87.9018, 0.0, 87.4182, 15.1921]
	c1 = [cb['radius']+414.0, 0.0006189, 51.6393, 0.0, 234.1955, 105.6372]
	c2 = [cb['radius']+35000.0, 0.0006189, 51.6393, 0.0, 234.1955, 105.6372]

	op0 = OP(c0,tspan,dt,coes=True)
	op1 = OP(c1,tspan,dt,coes=True)
	op2 = OP(c2,tspan,dt,coes=True)

	op0.propagate_orbit()
	op1.propagate_orbit()
	op2.propagate_orbit()

	t.plot_n_orbits([op0.rs, op1.rs, op2.rs], labels=['OneWeb', 'Iss', 'GEO'], show_plot=True)


	