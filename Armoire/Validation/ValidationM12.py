import numpy as np

from sys import path
path.append('D:\\Users\\s606482\\Desktop\\pyOrbit-master\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import planetary_data as pd
import tools as t



dt = 3600
cb = pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2']=True

    # MISSION 2
    
    perts['thrust']=0.045
    perts['isp']=1400
    mass0 = 1430
    tspan = dt*90*24
    # a,e,i,ta,aop,raan, date
    coes = np.array([cb['radius']+700, 0.005, 90, -106, 170, 1.6, [0,0,0]])
    op = OP(coes, tspan, dt, mass0, coes=True, deg=True, perts=perts)
    op.calculate_coes()
    op.plot_coes(days=True, show_plot=True)

    for i,coes in enumerate(op.coes):
        if coes[0]>1500+cb['radius']:
            print(op.ts[i]/(3600*24*360))
            break
