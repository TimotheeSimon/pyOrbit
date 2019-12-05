import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import math as ma

import planetary_data as pd
import tools as t

def null_perts():
    return {
        'J2':False,
        'aero':False,
        'moon_grav':False,
        'solar_grav':False,
        'thrust':0,
        'isp':0,
        'mdot':0,
        'direction':1,
        'a_maneuvre':True,
        'ecc_maneuvre':False,
        'inc_maneuvre':False,
        'raan_maneuvre':False,
        'aop_maneuvre':False,
    }

class OrbitPropagator:

    def __init__(self, state0, tspan, dt,
                 mass0=0,
                 t0=0,
                 coes=False,
                 deg=True,
                 cb=pd.earth,
                 perts=null_perts(),
                 propagator='lsoda'):
        
        if coes:
            self.r0, self.v0, self.date = t.coes2rv(state0, mu=cb['mu'], deg=deg)
                    
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]
                    
        self.tspan = tspan
        self.dt = dt
        self.cb = cb
        self.mass0 = mass0
        self.propagator = propagator

        self.n_steps = int(np.ceil(self.tspan/self.dt))

        self.ts = np.zeros((self.n_steps,1))
        self.ys = np.zeros((self.n_steps, 7))
        self.y0 = self.r0.tolist()+self.v0.tolist()+[self.mass0]
        self.ts[0] = t0
        self.ys[0] = self.y0
        self.step = 1

        self.solver = ode(self.diffy_q)
        self.solver.set_integrator(propagator)
        self.solver.set_initial_value(self.y0, t0)

        self.perts = perts

        self.propagate_orbit()

    def propagate_orbit(self):

        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step]=self.solver.t
            self.ys[self.step]=self.solver.y
            self.step+=1

        self.rs = self.ys[:,:3]
        self.vs = self.ys[:,3:6]
        self.masses = self.ys[:,-1]

    def diffy_q(self, time, y):
        rx,ry,rz,vx,vy,vz,mass = y
        
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])
        norm_r = np.linalg.norm(r)
        norm_v = np.linalg.norm(v)

        # Two bodies acceleration
        a = -r*self.cb['mu']/norm_r**3

        # J2 Perturbation
        if self.perts['J2']:
            h = np.cross(r,v)
            norm_h = np.linalg.norm(h)
            sma, e_norm, i, ta, aop, raan = t.rv2coes(r, v, mu=self.cb['mu'])
            
            A = - 1.5 * self.cb['J2'] * self.cb['mu'] * self.cb['radius']**2/norm_r**4
    
            pr = (1 - 3*ma.sin(i)**2*ma.sin(aop+ta)**2) * A * r/norm_r
            pv = ma.sin(i)**2*ma.sin(2*(aop+ta)) * A * v/norm_v
            ph = ma.sin(2*i)*ma.sin(aop+ta) * A * h/norm_h

            a_j2 = pr + pv + ph
    
            a += a_j2

        # Thrust

        mdot = - self.perts['mdot']

        if self.perts['thrust']:
            if not self.perts['J2']:
                h = np.cross(r,v)
                norm_h = np.linalg.norm(h)
                sma, e_norm, i, ta, aop, raan = t.rv2coes(r, v, mu=self.cb['mu'])

            c = np.cross(h,r)
            norm_c = np.linalg.norm(c)
            alpha, beta = 0, 0
            
            if self.perts['a_maneuvre']:
                alpha = ma.atan(e_norm*ma.sin(ta)/(1+e_norm*ma.sin(ta)))
                beta = 0

            if self.perts['ecc_maneuvre']:
                cosE = (e_norm+ma.cos(ta))/(1+e_norm*ma.cos(ta))
                alpha = ma.atan(ma.sin(ta)/(ma.cos(ta)+cosE))
                beta = 0

            if self.perts['inc_maneuvre']:
                alpha = 0
                beta = np.sign(ma.cos(aop+ta))*ma.pi/2

            if self.perts['raan_maneuvre']:
                alpha = 0
                beta = np.sign(ma.sin(aop+ta))*ma.pi/2

            if self.perts['aop_maneuvre']:
                alpha = ma.atan((1+e_norm*ma.cos(ta))/((2+e_norm*ma.cos(ta))*ma.tan(ta)))
                beta = (ma.atan(e_norm*ma.sin(aop+ta)/
                                (ma.tan(i)*(ma.sin(alpha-ta)*(1+e_norm*ma.cos(ta)))-ma.cos(alpha)*ma.sin(ta))))

            T = (self.perts['direction']*self.perts['thrust']*
            (ma.cos(beta)*ma.sin(alpha)*r/norm_r
             + ma.cos(beta)*ma.cos(alpha)*c/norm_c
             + ma.sin(beta)*h/norm_h)/1000)

            if not self.perts['mdot']:
                mdot = - self.perts['thrust']/self.perts['isp']/9.81

            a += T/mass
                         
        return [vx,vy,vz,a[0],a[1],a[2], mdot]

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

        axs[0,1].plot(ts,np.round(self.coes[:,1]),5)
        axs[0,1].set_title('Eccentricity vs. Time')
        axs[0,1].grid(True)
        

        axs[0,2].plot(ts,self.coes[:,4])
        axs[0,2].set_title('Argument of Periapse vs. Time')
        axs[0,2].grid(True)
        axs[0,2].set_ylabel('Angle (degrees)')

        axs[1,1].plot(ts,np.round(self.coes[:,2],5))
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

        if save_plot:
            plt.savefig(title+'.png', dpi=300)

    def plot_3d(self, cb=pd.earth,
                title='Figure',
                show_plot=True,
                save_plot=False,
                return_ax=False):
        
        fig =  plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111, projection='3d')
        #ax.set_aspect("equal")
        
        ax.plot(self.rs[:,0], self.rs[:,1], self.rs[:,2], color='xkcd:crimson', label="Computed trajectory")
        ax.plot([self.rs[0,0]], [self.rs[0,1]], [self.rs[0,2]], 'o', color='xkcd:green')


        _u, _v = np.mgrid[0:2*np.pi:25j, 0:np.pi:15j]
        _x = self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y = self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z = self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x, _y, _z, cmap='Blues')

        l = self.cb['radius']*1.5
        x,y,z=[[0,0,0], [0,0,0], [0,0,0]]
        u,v,w=[[l,0,0], [0,l,0], [0,0,l]]

        ax.quiver(x,y,z,u,v,w, color='k', arrow_length_ratio=0.1)

        max_val=np.max(np.abs(self.rs))
        ax.set_xlim(-max_val, max_val)
        ax.set_ylim(-max_val, max_val)
        ax.set_zlim(-max_val, max_val)

        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')

        plt.legend()
        plt.title(title)
        
        if show_plot:
            plt.show()
            
        if save_plot:
            plt.savefig(title+'.png', dpi=300)
            
        if return_ax:
            return ax

    def plot_masses(self, hours=False,
                    title='Masses',
                    days=False,
                    show_plot=True,
                    save_plot=False):
        
        print('Plotting Masses..')
        print('Initial mass :', self.masses[0],' kg')
        print('Final mass :', round(self.masses[-1],1),' kg')
        print('Consumed mass :', round(self.masses[0]-self.masses[-1],1),' kg')

        plt.figure(figsize=(10,6))

        if hours:
            ts = self.ts/3600.0
            xlabel = 'Time Elapsed (hours)'
        elif days:
            ts = self.ts/3600.0/24.0
            xlabel = 'Time Elapsed (days)'
        else:
            ts = self.ts
            xlabel = 'Time Elapsed (seconds)'

        plt.plot(ts, self.masses)
        plt.xlabel(xlabel)
        plt.ylabel('Mass (kg)')
        plt.title('Mass vs. Time')

        if show_plot:
            plt.show()

        if save_plot:
            plt.savefig(title+'.png', dpi=300)
            
        

        

        









