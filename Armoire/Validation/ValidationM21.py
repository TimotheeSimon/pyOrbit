import numpy as np

from sys import path
path.append('D:\\Users\\s606482\\Desktop\\pyOrbit-master\\pyOrbit')
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import planetary_data as pd
import tools as t



dt = 2*3600
cb = pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2']=True

    # MISSION 3
    
    perts['thrust']=0.09
    perts['isp']=1660
    mass0 = 1900
    tspan = dt*12*360*2
    # a,e,i,ta,aop,raan, date
    coes = np.array([cb['radius']+700, 0.005, 45, -106, 170, 1.6, [0,0,0]])
    op = OP(coes, tspan, dt, mass0, coes=True, deg=True, perts=perts)
    op.calculate_coes()
    op.plot_coes(days=True, show_plot=True)

    for i,coes in enumerate(op.coes):
        if coes[0]>8000+cb['radius']:
            print(op.ts[i]/(3600*24*360))
            break
