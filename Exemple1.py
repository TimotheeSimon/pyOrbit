import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

from sys import path
path.append('C:\\Users\\timsi\\Documents\\Visual Studio Code\\Orbit Visualisation\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
import planetary_data as pd

if __name__ == '__main__':
	cb = pd.earth
	altitude = cb['radius']+35000
	velocity = np.sqrt(cb['mu']/altitude)*0.7


	r0 = [altitude, 0, 0]
	v0 = [0, velocity, 0]

	dt = 60
	tspan = 24*60*dt

	op = OP(r0, v0, tspan, dt)
	op.propagate_orbit()
	op.plot_3d(show_plot=True)
	