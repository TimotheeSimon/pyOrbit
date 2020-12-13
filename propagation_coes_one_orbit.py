import numpy as np
import pandas as p
import matplotlib.pyplot as plt
import datetime



from sys import path
path.append('.\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import perturbations
import planetary_data as pd
import tools as t

            
    

dt = 60
tspan = 24*3600

if __name__ == '__main__':
    perts = perturbations()
    perts['J2']=True

    header = ['date', 'x', 'y', 'z', 'Vx', 'Vy', 'Vz', 'Ax', 'Ay', 'Az', 'a', 'e', 'i', 'aop', 'raan', 'M']
    data = p.read_csv('data_coes_one_orbit.csv', sep=';', names=header)
    data.date = p.to_datetime(data.date)


    Y0 = np.array([data.a[0], data.e[0], data.i[0], data.M[0], data.aop[0], data.raan[0], [0, 0, 0]])


    Hispasat_sim = OP(Y0, tspan, dt, deg=True, coes=True, perts=perts)
    Hispasat_sim.propagate_orbit()
    Hispasat_sim.calculate_coes()
    
    coes = Hispasat_sim.coes
    
    T_sim = [p.to_datetime(t, unit='s', origin=data.date[0]) for t in Hispasat_sim.ts]

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
    
    # a, ecc, inc, ta, argp, raan
    
    axs[0,0].plot(T_sim, coes[:,0])
    axs[0,0].plot(data.date, data.a, '.', markersize=0.5)
    axs[0,0].set_title('Demi grand axe (km)')
    axs[0,0].set_xticks([])
    
    axs[0,1].plot(T_sim, coes[:,1])
    axs[0,1].plot(data.date, data.e, '.', markersize=0.5)
    axs[0,1].set_title('Eccentricité')
    axs[0,1].set_xticks([])
    
    axs[0,2].plot(T_sim, coes[:,2])
    axs[0,2].plot(data.date, data.i, '.', markersize=0.5)
    axs[0,2].set_title('Inclinaison (deg)')
    axs[0,2].set_xticks([])

    axs[1,0].plot(T_sim, coes[:,4])
    axs[1,0].plot(data.date, data.aop, '.', markersize=0.5)
    axs[1,0].set_title('Argument du périgée (deg)')

    axs[1,1].plot(T_sim, coes[:,5])
    axs[1,1].plot(data.date, data.raan, '.', markersize=0.5)
    axs[1,1].set_title('Longitude du noeud ascendant (deg)')

    axs[1,2].plot(T_sim, coes[:,3], label='Simulé')
    axs[1,2].plot(data.date, data.M, '-', markersize=0.5, label='Mesuré')
    axs[1,2].set_title('Anomalie moyenne (deg)')

    for ax in fig.axes:
        plt.sca(ax)
        plt.xticks(rotation=90)

        
    plt.legend()
    plt.subplots_adjust(hspace=0.25)
    plt.show()
    
