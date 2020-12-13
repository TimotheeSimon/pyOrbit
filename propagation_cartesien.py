import numpy as np
from math import sqrt
from PyAstronomy import pyasl
import pandas as p
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
import datetime



from sys import path
path.append('.\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import perturbations
import planetary_data as pd
import tools as t
    

dt = 60
tspan = 8*24*3600

if __name__ == '__main__':
    perts = perturbations()
    perts['J2'] = True
    header = ['date', 'x', 'y', 'z', 'Vx', 'Vy', 'Vz', 'Ax', 'Ay', 'Az']
    data = p.read_csv('data_cartesien.csv', sep=';', names=header)
    data.date = p.to_datetime(data.date)

    Y0 = np.array([data.x[0], data.y[0], data.z[0], data.Vx[0], data.Vy[0], data.Vz[0]])


    Hispasat_sim = OP(Y0, tspan, dt, perts=perts)
    Hispasat_sim.propagate_orbit()

    pos = Hispasat_sim.rs
    vit = Hispasat_sim.vs
    T_sim = [p.to_datetime(t, unit='s', origin=data.date[0]) for t in Hispasat_sim.ts]

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
    
    
    axs[0,0].plot(T_sim, pos[:,0])
    axs[0,0].plot(data.date, data.x, '.', markersize=0.5)
    axs[0,0].set_title('X')
    axs[0,0].set_ylabel('Position')
    axs[0,0].set_xticks([])
    
    axs[0,1].plot(T_sim, pos[:,1])
    axs[0,1].plot(data.date, data.y, '.', markersize=0.5)
    axs[0,1].set_title('Y')
    axs[0,1].set_xticks([])
    
    axs[0,2].plot(T_sim, pos[:,2])
    axs[0,2].plot(data.date, data.z, '.', markersize=0.5)
    axs[0,2].set_title('Z')
    axs[0,2].set_xticks([])

    axs[1,0].plot(T_sim, vit[:,0])
    axs[1,0].plot(data.date, data.Vx, '.', markersize=0.5)
    axs[1,0].set_title('Vx')

    axs[1,1].plot(T_sim, vit[:,1])
    axs[1,1].plot(data.date, data.Vy, '.', markersize=0.5)
    axs[1,1].set_title('Vy')

    axs[1,2].plot(T_sim, vit[:,2], label='Simulé')
    axs[1,2].plot(data.date, data.Vz, '.', markersize=0.5, label='Mesuré')
    axs[1,2].set_title('Vz')
    

    for ax in fig.axes:
        plt.sca(ax)
        plt.xticks(rotation=90)
    
    plt.legend()
    plt.subplots_adjust(hspace=0.25)
    plt.show()
    
