import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import math as ma
import vpython as v

import planetary_data as pd
import tools as t

def perturbations():

    """ Défini les perturbations et tous les paramètre qui leur sont relatifs.
    à noter : cf. examples pour comprendre le fonctionnement

    outputs:
    -- dictionnaire python contenant tous les paramètres relatifs aux perturbations
    """

    return {
        'J2':False,
        'aero':False,
        'moon_grav':False,
        'solar_grav':False,
        'thrust':0,
        'isp':0,
        'mdot':0,
        'direction':1,
        'a_maneuvre':False,
        'ecc_maneuvre':False,
        'inc_maneuvre':False,
        'raan_maneuvre':False,
        'aop_maneuvre':False,
        }


class OrbitPropagator:
    """ Cette classe défini l'object permettant la propagation d'une orbite pour un état initial donné
    """

    def __init__(self, state0, tspan, dt,
                 mass0=0,
                 t0=0,
                 coes=False,
                 deg=True,
                 cb=pd.earth,
                 perts=perturbations(),
                 propagator='vode'):


        """ Constructeur de classe
        inputs:
        state0 - vecteur d'état initial - array 1x6 [km, km, km, km/s, km/s, km/s] ou array 1x7 [km, 1, rad, rad, rad, rad, [année, mois, jours, heures]]
        mass0 - masse initiale - float [kg] (par défaut : 0)
        t0 - instant initial - float [s] (par défaut : 0)
        coes - Indique si le vecteur d'état initial contient les elements Keplerien les vecteurs position/vitesse - boolean (par défaut : False)
        deg - Dans le cas ou les elements Keplerien sont utilisés comme entrée, indique s'ils sont exprimés en degrés ou en radians - boolean (par défaut : True)
        cb - dictionnaire python définissant le corps central d'une orbite - (par défaut : la terre)
        perts - dictionnaire python définissant les perturbations - (par défaut : cf. fonction perturbations())
        propagator - solver utilisé lors de la propagation de l'orbite - string (par défaut : 'dop853' !!!!!!! plus précis que 'lsoda' mais plus lent)
        à noter : "dop853" est en fait un RK 8

        """
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
    
        """Fonction appelée directement dans le constructeur
        propage l'orbite dans le temps en enregistrant les etats a chaque itération
        """

        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step]=self.solver.t
            self.ys[self.step]=self.solver.y
            self.step+=1

        self.rs = self.ys[:,:3]
        self.vs = self.ys[:,3:6]
        self.masses = self.ys[:,-1]

    def diffy_q(self, time, y):

        """calcul dX/dt avec X le vecteur d'état

        inputs:
        time -- temps courant - float [s]
        y -- vecteur d'état courant - array of float [km, km, km, km/s, km/s, km/s, kg]

        outputs:
        [vx,vy,vz,a[0],a[1],a[2], mdot] -- dérivée du vecteur d'état - array of float [km/s, km/s, km/s, km/s², km/s², km/s², kg/s]
        """
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

            try:
                a += a_j2
            except:
                print("Warning: mass must be set in constructor")

        # Thrust

        # Dans l'état des choses, les résultats sont satisfaisant pour des manoeuvres ayant pour but la modification d'un seul paramètre à la fois.
        # l'orientation du vecteur poussée est calculée selon la méthode expliquée dans "Ressources/Comparison_of_Low_Thrust_Control_Laws.pdf" et "Ressources/Low_Thrust_Maneuvers.pdf"
        # La modification de l'eccentricité ne semble pas fonctionner (cf. "Exploitation/Examples/thrust_ecc.py")


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

        """ calcul tous les elements Keplerien après avoir propagé l'orbite sur les vecteurs position et vitesse
        les elements Kepleriens calculés sont stockés dans l'attribut "coes"  de l'objet OrbitPropagator sous forme d'un array nx7.

        inputs:
        deg -- impose les unité des angles des paramètres Képlériens - boolean (par defaut True)
        """

        print('Calcutating classical orbital elements ..')

        self.coes = np.zeros((self.n_steps,6))

        for n in range(self.step):
            self.coes[n,:] = t.rv2coes(self.rs[n,:], self.vs[n,:], mu=self.cb['mu'], deg=deg)

    def plot_coes(self, hours=False, days=False, show_plot=False, save_plot=False, title='COEs'):

        """ Trace l'évolution des elements Kepleriens en fonction du temps

        inputs:
        hours -- impose l'heure comme unité temporelle de tracé - boolean (par defaut False)
        days -- impose le jour comme unité temporelle de tracé - boolean (par defaut False)
        show_plot -- impose ou non le tracé à la fin de l'exécution de la fonction - boolean (par defaut False)
        save_plot -- impose ou non la sauvegrade du graphe à la fin de l'exécution de la fonction - boolean (par defaut False)
        title -- titre du graphe - string (par defaut 'COEs')
        à noter : par défaut, l'unité temporelle du tracé est la seconde.
        """

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


        axs[0,2].plot(ts,self.coes[:,4])
        axs[0,2].set_title('Argument of Periapse vs. Time')
        axs[0,2].grid(True)
        axs[0,2].set_ylabel('Angle (degrees)')

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

        if save_plot:
            plt.savefig(title+'.png', dpi=300)

    def plot_3d(self, cb=pd.earth,
                title='Figure',
                show_plot=True,
                save_plot=False,
                return_ax=False):

        """ Trace l'orbite dans l'espace

        inputs:
        cb - dictionnaire python définissant le corps central d'une orbite - (par défaut : la terre)
        title -- titre du graphe - string (par defaut 'Figure')
        show_plot -- impose ou non le tracé à la fin de l'exécution de la fonction - boolean (par defaut False)
        save_plot -- impose ou non la sauvegrade du graphe à la fin de l'exécution de la fonction - boolean (par defaut False)
        return_ax -- impose lerenvoie de l'objet "axe" à la fin de l'execution de la fonction - boolean (par défaut : False)

        outputs:
        (facultatif) ax -- objet matplotlib Axe
        """

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

    def plot_3d_animation(self):

        """ affiche dans le navigateur une animation dynamique de l'orbite
        """

        scene = v.canvas(width=1500, height=700, center=v.vector(0,5,0), background=v.color.white)
        scene.lights = []
        scene.ambient = v.color.gray(0.8)
        v.sphere(pos=v.vector(0,0,0), radius=self.cb['radius'], texture=v.textures.earth)
        r0 = self.rs[0]
        spacecraft = v.sphere(pos=v.vector(r0[0],r0[1],r0[2]), radius=250, color=v.vector(0, 21/25, 0), trail_color=v.vector(21/25, 0, 0), make_trail=True, trail_radius=20)

        for r in self.rs:
            v.rate(1000)
            spacecraft.pos = v.vector(r[0],r[1],r[2])





    def plot_masses(self, hours=False,
                    days=False,
                    show_plot=True,
                    save_plot=False,
                    title='Masses'):

        """ Trace l'évolution de la masse du satellite en fonction du temps

        inputs:
        hours -- impose l'heure comme unité temporelle de tracé - boolean (par defaut False)
        days -- impose le jour comme unité temporelle de tracé - boolean (par defaut False)
        show_plot -- impose ou non le tracé à la fin de l'exécution de la fonction - boolean (par defaut False)
        save_plot -- impose ou non la sauvegrade du graphe à la fin de l'exécution de la fonction - boolean (par defaut False)
        title -- titre du graphe - string (par defaut 'Masses')
        à noter : par défaut, l'unité temporelle du tracé est la seconde.
        """

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















