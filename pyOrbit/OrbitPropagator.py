import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

import planetary_data as pd
import tools as t

def null_perts():
	return {
			'J2':False,
			'aero':False,
			'moon_grav':False,
			'solar_grav':False
	}

class OrbitPropagator:

	def __init__(self, state0, tspan, dt, coes=False,cb=pd.earth, deg=True, perts=null_perts()):
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

		self.perts = perts

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

		# Two bodies acceleration
		a = -r*self.cb['mu']/norm_r**3

		# J2 Perturbation
		if self.perts['J2']:
			z2 = r[2]**2
			r2 = norm_r**2
			tx = r[0]/norm_r*(5*z2/r2-1)
			ty = r[1]/norm_r*(5*z2/r2-1)
			tz = r[2]/norm_r*(5*z2/r2-3)

			a_j2 = 1.5 * self.cb['J2'] * self.cb['mu'] * self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])
			a+=a_j2


		return[vx,vy,vz,a[0],a[1],a[2]]

	def calculate_coes(self, deg=True):
		print('Calcutating classical orbital elements ..')

		self.coes = np.zeros((self.n_steps,6))
		
		for n in range(self.step):
			self.coes[n,:] = t.rv2coes(self.rs[n,:], self.vs[n,:], mu=self.cb['mu'], deg=deg)
	
	def plot_coes(self, hours=False, days=False, show_plot=False, save_plot=False, title='COEs'):
		print('Plotting COEs..')

		fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15,9))

		fig.suptitle(title, fontsize=20)

		if hours:
			ts = self.ts/3600.0
			xlabel = 'Time Elapsed (hours)'
		elif days:
			ts = self.ts/3600.0/24.0
			xlabel = 'Time Elapsed (days)'
		else:
			ts = self.ts
			xlabel = 'Time Elapsed (seconds)'
		
		axs[0,0].plot(ts,self.coes[:,3])
		axs[0,0].set_title('True Anomaly vs. Time')
		axs[0,0].grid(True)
		axs[0,0].set_ylabel('Angle (degrees)')

		axs[1,0].plot(ts,self.coes[:,0])
		axs[1,0].set_title('Semi-Major Axis vs. Time')
		axs[1,0].grid(True)
		axs[1,0].set_ylabel('Semi-Major Axis (km)')
		axs[1,0].set_xlabel(xlabel)

		axs[0,1].plot(ts,self.coes[:,1])
		axs[0,1].set_title('Eccentricity vs. Time')
		axs[0,1].grid(True)
		axs[0,1].set_ylabel('Angle (degrees)')

		axs[0,2].plot(ts,self.coes[:,4])
		axs[0,2].set_title('Argument of Periapse vs. Time')
		axs[0,2].grid(True)

		axs[1,1].plot(ts,self.coes[:,2])
		axs[1,1].set_title('Inclinaison vs. Time')
		axs[1,1].grid(True)
		axs[1,1].set_ylabel('Angle (degrees)')
		axs[1,1].set_xlabel(xlabel)

		axs[1,2].plot(ts,self.coes[:,5])
		axs[1,2].set_title('RAAN vs. Time')
		axs[1,2].grid(True)
		axs[1,2].set_xlabel(xlabel)

		plt.subplots_adjust(hspace=0.3, wspace=0.4)

		if show_plot:
			plt.show()

		if show_plot:
			plt.savefig(title+'.png', dpi=300)










