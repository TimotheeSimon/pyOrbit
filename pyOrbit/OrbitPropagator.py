import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

import planetary_data as pd
import tools as t

class OrbitPropagator:

	def __init__(self, state0, tspan, dt, coes=False,cb=pd.earth, deg=True):
		if coes:
			self.r0, self.v0, self.date = t.coes2rv(state0, mu=cb['mu'], deg=deg)
		else:
			self.r0 = state0[:3]
			self.v0 = state0[3:]
			
		self.tspan = tspan
		self.dt = dt
		self.cb = cb

		self.n_steps = int(np.ceil(self.tspan/self.dt))

		self.ts = np.zeros((self.n_steps,1))
		self.ys = np.zeros((self.n_steps, 6))
		self.y0 = self.r0.tolist()+self.v0.tolist()
		self.ts[0] = 0
		self.ys[0] = self.y0
		self.step = 1

		self.solver = ode(self.diffy_q)
		self.solver.set_integrator('lsoda')
		self.solver.set_initial_value(self.y0)

		self.propagate_orbit()
	
	def propagate_orbit(self):

		while self.solver.successful() and self.step<self.n_steps:
			self.solver.integrate(self.solver.t+self.dt)
			self.ts[self.step]=self.solver.t
			self.ys[self.step]=self.solver.y
			self.step+=1

		self.rs = self.ys[:,:3]
		self.vs = self.ys[:,3:]

	def diffy_q(self, t, y):
		rx,ry,rz,vx,vy,vz = y
		r = np.array([rx,ry,rz])
		v = np.array([vx,vy,vz])
		norm_r = np.linalg.norm(r)
		ax,ay,az = -r*self.cb['mu']/norm_r**3

		return[vx,vy,vz,ax,ay,az]

	

	








