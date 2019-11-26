import numpy as np
import math as ma
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import planetary_data as pd

d2r = np.pi/180


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
    a,e,i,ta,aop,raan = coes
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

    return r, v

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