import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D



from sys import path
path.append('.\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import planetary_data as pd
import tools as t




cb = pd.earth

if __name__ == '__main__':

    fichier = open("manoeuvre.txt")

    a = []
    ecc = []
    inc = []
    raan = []
    argp = []

    for line in fichier.readlines():
        mess = line.strip().split('\t')
        a.append(float(mess[0]))
        ecc.append(float(mess[1]))
        inc.append(float(mess[2]))
        raan.append(float(mess[3]))
        argp.append(float(mess[4]))
    fichier.close()

    a = np.array(a)
    ecc = np.array(ecc)
    inc = np.array(inc)
    raan = np.array(raan)*180/np.pi
    argp = np.array(argp)

    
    perts = null_perts()
    perts['J2']=False
    perts['thrust']=0
    perts['isp']=1450
    perts['direction']=-1
    ## Orbite passive initiale
    
    
    T = []
    t0 = 0
    masses = []
    mass0 = 2000
    dt = 60
    tspan = 3600*24*19
    dt = 3600
    i0 = 120
    
    coes = np.array([a[i0], ecc[i0], inc[i0], 0, argp[i0], raan[i0], [0,0,0]])
    op = OP(coes, tspan, dt, mass0, t0, coes=True, deg=True, perts=perts)
    op.calculate_coes()
    
    Y = op.coes
    T = [float(t) for t in op.ts]

    ## Poussée pour réduire l'altitude

    coes0 = op.coes[-1,:]
    t0 = T[-1]
    mass0 = op.masses[-1]
    dt = 60
    tspan = 3600*24

    perts['thrust']=0.085*2
    
    coes = np.array([coes0[0], coes0[1], coes0[2], coes0[3], coes0[4], coes0[5], [0,0,0]])
    op = OP(coes, tspan, dt, mass0, t0, coes=True, deg=True, perts=perts)
    op.calculate_coes()

    Y = np.append(Y,op.coes, axis=0)
    T = np.append(T,[float(t) for t in op.ts])

    ## Orbite passive transitoire 1
    dt = 120
    
    coes0 = op.coes[-1,:]
    t0 = T[-1]
    mass0 = op.masses[-1]
    tspan = 3600*24*37 - 3600
    
    perts['thrust']=0
    perts['direction']=1

    a0 = op.coes[0,0]
    

    coes = np.array([coes0[0], coes0[1], coes0[2], coes0[3], coes0[4], coes0[5], [0,0,0]])
    op = OP(coes, tspan, dt, mass0, t0, coes=True, deg=True, perts=perts)
    op.calculate_coes()

    Y = np.append(Y,op.coes, axis=0)
    T = np.append(T,[float(t) for t in op.ts])
    ## Poussée pour augmenter l'altitude
    dt = 60*15
    
    coes0 = op.coes[-1,:]
    t0 = T[-1]
    mass0 = op.masses[-1]
    tspan = 3600*24*2

    perts['thrust']=0.018
    perts['direction']=1
    
    coes = np.array([coes0[0], coes0[1], coes0[2], coes0[3], coes0[4], coes0[5], [0,0,0]])
    op = OP(coes, tspan, dt, mass0, t0, coes=True, deg=True, perts=perts)
    op.calculate_coes()

    Y = np.append(Y,op.coes, axis=0)
    T = np.append(T,[float(t) for t in op.ts])

    ## Orbite passive transitoire 2
    dt = 3600
    
    coes0 = op.coes[-1,:]
    t0 = T[-1]
    mass0 = op.masses[-1]
    tspan = 3600*24*5
    
    perts['thrust']=0

    coes = np.array([coes0[0], coes0[1], coes0[2], coes0[3], coes0[4], coes0[5], [0,0,0]])
    op = OP(coes, tspan, dt, mass0, t0, coes=True, deg=True, perts=perts)
    op.calculate_coes()

    Y = np.append(Y,op.coes, axis=0)
    T = np.append(T,[float(t) for t in op.ts])

    ## Poussée pour augmenter l'altitude 2
    dt = 60*10
    
    coes0 = op.coes[-1,:]
    t0 = T[-1]
    mass0 = op.masses[-1]
    tspan = 3600*(26*2 - 5) 

    perts['thrust']=0.069
    perts['direction']=1
    
    coes = np.array([coes0[0], coes0[1], coes0[2], coes0[3], coes0[4], coes0[5], [0,0,0]])
    op = OP(coes, tspan, dt, mass0, t0, coes=True, deg=True, perts=perts)
    op.calculate_coes()

    Y = np.append(Y,op.coes, axis=0)
    T = np.append(T,[float(t) for t in op.ts])

    ## Orbite passive transitoire 3
    dt = 3600
    
    coes0 = op.coes[-1,:]
    t0 = T[-1]
    mass0 = op.masses[-1]
    tspan = 3600*24*21
    perts['thrust']=0

    coes = np.array([coes0[0], coes0[1], coes0[2], coes0[3], coes0[4], coes0[5], [0,0,0]])
    op = OP(coes, tspan, dt, mass0, t0, coes=True, deg=True, perts=perts)
    op.calculate_coes()

    Y = np.append(Y,op.coes, axis=0)
    T = np.append(T,[float(t) for t in op.ts])

    


    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(14,10))
    


    
    axs[0].plot(np.array(T)/3600/24,Y[:,0])
    axs[0].plot(a[i0:])
    axs[0].plot(a[i0:], '.')
    axs[0].set_title('Demi grand axe vs. Temps')
    axs[0].set_ylabel('Demi grand axe (km)')
 

    axs[1].plot(np.array(T)/3600/24,Y[:,1])
    axs[1].plot(ecc[i0:])
    axs[1].plot(ecc[i0:], '.')
    axs[1].set_title('Eccentricité vs. Temps')
    axs[1].set_ylabel('Eccentricité')


    axs[2].plot(np.array(T)/3600/24,Y[:,4])
    axs[2].plot(argp[i0:])
    axs[2].plot(argp[i0:], '.')
    axs[2].set_title('Argument du périgée vs. Temps')
    axs[2].set_xlabel('Temps (jours)')
    axs[2].set_ylabel('Argument du périgée (deg)')

    plt.subplots_adjust(hspace=0.25)
    plt.show()


    
