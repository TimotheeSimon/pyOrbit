import numpy as np
import math as ma
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


import planetary_data as pd

d2r = np.pi/180

def plot_3d(rs, cb=pd.earth, title='Figure', show_plot=False, save_plot=False):
		fig =  plt.figure(figsize=(10,8))
		ax = fig.add_subplot(111, projection='3d')
		#ax.set_aspect("equal")
		
		ax.plot(rs[:,0], rs[:,1], rs[:,2], color='xkcd:crimson', label="Trajectory")
		ax.plot([rs[0,0]], [rs[0,1]], [rs[0,2]], 'o', color='xkcd:crimson')


		_u, _v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
		_x = cb['radius']*np.cos(_u)*np.sin(_v)
		_y = cb['radius']*np.sin(_u)*np.sin(_v)
		_z = cb['radius']*np.cos(_v)
		ax.plot_surface(_x, _y, _z, cmap='Blues')

		l=cb['radius']*1.5
		x,y,z=[[0,0,0], [0,0,0], [0,0,0]]
		u,v,w=[[l,0,0], [0,l,0], [0,0,l]]

		ax.quiver(x,y,z,u,v,w, color='k', arrow_length_ratio=0.1)

		max_val=np.max(np.abs(rs))
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

def plot_n_orbits(rs, labels, cb=pd.earth, title='Figure', show_plot=False, save_plot=False):
    fig =  plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection='3d')
    
    
    n=0
    for r in rs:
        ax.plot(r[:,0], r[:,1], r[:,2], label=labels[n])
        ax.plot([r[0,0]], [r[0,1]], [r[0,2]], 'o',)
        n+=1


    _u, _v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
    _x = cb['radius']*np.cos(_u)*np.sin(_v)
    _y = cb['radius']*np.sin(_u)*np.sin(_v)
    _z = cb['radius']*np.cos(_v)
    ax.plot_surface(_x, _y, _z, cmap='Blues')

    l=cb['radius']*1.5
    x,y,z=[[0,0,0], [0,0,0], [0,0,0]]
    u,v,w=[[l,0,0], [0,l,0], [0,0,l]]

    ax.quiver(x,y,z,u,v,w, color='k', arrow_length_ratio=0.1)

    max_val=np.max(np.abs(rs))
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

def coes2rv(coes, mu=pd.earth['mu'], deg=True):
    a,e,i,ta,aop,raan, date = coes
    if deg:
        i*=d2r
        ta*=d2r
        aop*=d2r
        raan*=d2r
    
    E = ecc_anomay([ta,e],'tae')
    # Ne serait-ce pas un '+', ci dessous ?
    r_norm = a*(1-e**2)/(1+e*np.cos(ta))

    # Périfocal
    r_perif = r_norm*np.array([ma.cos(ta),ma.sin(ta),0])
    v_perif = ma.sqrt(mu*a)/r_norm*np.array([-ma.sin(E),ma.cos(E)*ma.sqrt(1-e**2),0])

    perif2eci = np.transpose(eci2perif(raan,aop,i))

    r = np.dot(perif2eci, r_perif)
    v = np.dot(perif2eci, v_perif)

    return r, v, date

def ecc_anomay(arr, method, tol=1e-8, max_step=200):
    if method=='newton':
        Me,e=arr
        if Me<np.pi/2.0: E0=Me+e/2
        else: E0=Me-e
        for n in range(max_step):
            ratio=(E0-e*np.sin(E0)-Me)/(1-e*np.cos(E0))
            if abs(ratio)<tol:
                if n==0: return E0
                else: return E1
            E1 = E0-ratio
            E0=E1
        return False
    elif method=='tae':
        ta,e=arr
        return 2*ma.atan(ma.sqrt((1-e)/(1+e))*ma.tan(ta/2.0))
    else:
        print('Invalid method for eccentric anomaly')    

def eci2perif(raan,aop,i):


    row0 = [ma.cos(raan)*ma.cos(aop)-ma.sin(raan)*ma.sin(aop)*ma.cos(i), 
            ma.sin(raan)*ma.cos(aop)+ma.cos(raan)*ma.sin(aop)*ma.cos(i), 
            ma.sin(i)*ma.sin(aop)]
    row1 = [-ma.cos(raan)*ma.sin(aop)-ma.sin(raan)*ma.cos(aop)*ma.cos(i), 
            -ma.sin(raan)*ma.sin(aop)+ma.cos(raan)*ma.cos(aop)*ma.cos(i), 
            ma.sin(i)*ma.cos(aop)]    
    row2 = [ma.sin(raan)*ma.sin(i), 
            -ma.cos(raan)*ma.sin(i), 
            ma.cos(i)]    
    return np.array([row0, row1, row2])

def tle2coes(tle_filename, mu=pd.earth['mu']):
    with open(tle_filename, 'r') as f:
        lines = f.readlines()

    line0 = lines[0].strip()
    line1 = lines[1].strip().split()
    line2 = lines[2].strip().split()

    epoch = line1[3]
    year, month, day, hour = calc_epoch(epoch)

    i = float(line2[2])*d2r
    raan = float(line2[3])*d2r
    e_string = line2[4]
    e = float('0.'+e_string)
    aop = float(line2[5])*d2r
    Me = float(line2[6])*d2r
    mean_motion = float(line2[7]) # revs/day
    T = 1/mean_motion*24*3600 # seconds
    a = (T**2*mu/4.0/np.pi**2)**(1/3.0)

    E = ecc_anomay([Me,e], 'newton')
    ta = true_anomaly([E, e])
    r_mag = a*(e-np.cos(E))


    return a, e, i, ta, aop, raan, [year, month, day, hour]

def calc_epoch(epoch):
    year = int('20'+epoch[:2])

    epoch = epoch[2:].split('.')

    day_of_year = int(epoch[0])-1
    hour = float('0.'+epoch[1])*24.0
    date = datetime.date(year, 1, 1)+datetime.timedelta(day_of_year)

    month = float(date.month)
    day = float(date.day)

    return year, month, day, hour

def true_anomaly(arr):
    E, e = arr
    return 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2.0))

def tle2rv(tle_filename):
    return coes2rv(tle2coes(tle_filename))