import numpy as np
import math as ma
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
import planetary_data as pd

d2r = np.pi/180
r2d = 180/np.pi

def coes2rv(coes, mu=pd.earth['mu'], deg=True):

    """Converti les elements kepleriens en vecteurs position et vitesse
    à noter: la date n'est utile que pour les TLE. la fixer à [0, 0, 0, 0] ne change pas les resultats.

    inputs:
    coes -- vecteur contenant les elements kepleriens - array 1x7 [km, 1, rad, rad, rad, rad, [année, mois, jours, heures]]
    mu -- parametre gravitationnel (km3/s2) - float (par defaut 398695.9848)
    deg -- impose les unité des angles en sorties - boolean (par defaut True)
        
    outputs:
    r -- vecteur position - array 1x3 [km]
    v -- vecteur vitesse - array 1x3 [km/s]
    date -- vecteur contenant la date  - array 1x4 [jours, heures, minutes, secondes]
    """

    a,e,i,ta,aop,raan, date = coes
    if deg:
        i*=d2r
        ta*=d2r
        aop*=d2r
        raan*=d2r
   
    E = ecc_anomay([ta,e],'tae')
    
    r_norm = a*(1-e**2)/(1+e*np.cos(ta))

    # Périfocal
    r_perif = r_norm*np.array([ma.cos(ta),ma.sin(ta),0])
    v_perif = ma.sqrt(mu*a)/r_norm*np.array([-ma.sin(E),ma.cos(E)*ma.sqrt(1-e**2),0])

    perif2eci = np.transpose(eci2perif(raan,aop,i))

    r = np.dot(perif2eci, r_perif)
    v = np.dot(perif2eci, v_perif)

    return r, v, date

def rv2coes(R, V, mu=pd.earth['mu'], deg=False, print_results=False, tol=1e-7):

    """Converti les vecteurs positions et vitesse en elements kepleriens
    à noter : hautement inspiré de l'algorithme 4.1 dans "Ressources/Orbital_mecanics_for_engineering_students.pdf"

    inputs:
    R -- vecteur position - array 1x3 [km]
    V -- vecteur vitesse(km/s) - array 1x3 [km]
    mu -- parametre gravitationnel - float [km3/s2] (par defaut 398695.9848) 
    deg -- impose les unité des angles en sorties - boolean (par defaut False)
    print_results -- impose l'affichage des resultats à la fin de la conversion - boolean (par defaut False)
    tol -- donne une tolérance en dessous de laquelle l'eccentricité est considéré nulle - float (par défaut 1e-7)
        
    outputs:
    [a, e, i, ta, aop, raan] -- vecteur contenant les elements kepleriens - array 1x6 [km/1/rad/rad/rad/rad]
    """
        
    r = norm(R)
    v = norm(V)
    rv = np.dot(R,V)/r
    
    H = np.cross(R,V)
    h = norm(H)
    
    if(h!=0):
        i = ma.acos(H[2]/h)
    else:
        i = 0
    
    N = np.cross([0,0,1], H)
    n = norm(N)
    
    if n != 0:
        raan = ma.acos(N[0]/n)
        if N[1]<0:
            raan = 2*np.pi - raan
    else:
        raan = 0
    
    E = 1/mu*((v**2-mu/r)*R - r*rv*V)
    e = norm(E)
    
    if n != 0:
        if e > tol:
            try:
                aop = ma.acos(np.dot(N,E)/(n*e))
            except:
                # renvoie 1 ou -1
                # Permet de s'assurer qu'il n'y a pas d'erreur d'arrondi numérique
                aop = np.dot(N,E)/abs(np.dot(N,E))
            if E[2]<0:
                aop = 2*np.pi - aop
        else:
            print('Verifier valeur AOP')
            aop = ma.pi/2
    else:
        aop = 0
    
    if e > tol:
        try:
            ta = ma.acos(np.dot(E,R)/(e*r))
        except:
            ta = np.dot(E,R)/abs(np.dot(E,R))
        if rv < 0:
            ta =  2*np.pi - ta
    else:
        cp = np.cross(N,R)
        if cp[2] >= 0:
            ta = ma.acos(np.dot(N,R)/(n*r))
        else:
            ta = 2*np.pi - ma.acos(np.dot(N,R)/(n*r))
    
    a = h**2/mu/(1-e**2)
        
    if print_results and deg:
        print('a : '+str(a))
        print('e : '+str(e))
        print('i : '+str(i*r2d))
        print('raan : '+str(raan*r2d))
        print('aop : '+str(aop*r2d))
        print('ta : '+str(ta*r2d))

    if print_results and not deg:
        print('a : '+str(a))
        print('e : '+str(e))
        print('i : '+str(i))
        print('raan : '+str(raan))
        print('aop : '+str(aop))
        print('ta : '+str(ta))
    
    if deg: return [a, e, i*r2d, ta*r2d, aop*r2d, raan*r2d]
    else: return [a, e, i, ta, aop, raan]

def ecc_anomay(arr, method, tol=1e-8, max_step=200):

    """Resoud l'equation de Kepler pour renvoyer l'anomalie excentrique
    à noter: Hautement inspiré de l'algorithm 3.3 présenté dans "Ressources/Orbital_mechanics_for_engineering_students.pdf"

    inputs:
    arr -- vecteur contenant l'anomalie moyenne et l'eccentricité  - array 1x2 [rad, 1]
    method -- string indiquant le nom de la méthode de calcul - 'newton' ou 'tae'
    tol -- tolérance de l'algorithme - float [1] (par défaut 1e-8)
    max_step -- nombre maximum d'itération de l'algorithme - float [1] (par défaut 200)
    
    outputs:
    E0 -- Anomalie excentrique - float [1]
    """

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

def true_anomaly(arr):

    """Calcul l'anomalie vraie à partir de l'anomalie excentrique

    inputs:
    arr -- vecteur contenant l'anomalie excentrique et l'eccentricité  - array 1x2 [rad, 1]

    outputs:
    ta -- anomalie vraie - float [rad]
    """
    
    E, e = arr
    return 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2.0))

def eci2perif(raan,aop,i):

    """Renvoie la matrice de passage du référentiel périfocal vers le référentiel ECI

    inputs:
    raan -- longitude du noeud ascendant - float [rad]
    aop -- argument du périgée - float [rad]
    i -- inclinaison - float [rad]
    
    outputs:
    M -- matrice de passage - array numpy 3x3 
    """
    
    row0 = [ma.cos(raan)*ma.cos(aop)-ma.sin(raan)*ma.sin(aop)*ma.cos(i), 
            ma.sin(raan)*ma.cos(aop)+ma.cos(raan)*ma.sin(aop)*ma.cos(i), 
            ma.sin(i)*ma.sin(aop)]
    row1 = [-ma.cos(raan)*ma.sin(aop)-ma.sin(raan)*ma.cos(aop)*ma.cos(i), 
            -ma.sin(raan)*ma.sin(aop)+ma.cos(raan)*ma.cos(aop)*ma.cos(i), 
            ma.sin(i)*ma.cos(aop)]    
    row2 = [ma.sin(raan)*ma.sin(i), 
            -ma.cos(raan)*ma.sin(i), 
            ma.cos(i)]    
    M = np.array([row0, row1, row2])
    return M
    
def tle2coes(tle_filename, mu=pd.earth['mu'], deg=False):

    """Extrait les elements kepleriens contenus dans un fichier TLE

    inputs:
    tle_filename -- chemin relatif et nom du fichier - string
    mu -- parametre gravitationnel (km3/s2) - float (par defaut 398695.9848)
    
    outputs:
    [a, e, i, ta, aop, raan, [year, month, day, hour]] -- vecteur contenant les elements kepleriens - array 1x7 [km, 1, rad, rad, rad, rad, [année, mois, jours, heures]]
    """
    
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

    E = ecc_anomay([Me,e], 'tae')
    ta = true_anomaly([E, e])
    r_mag = a*(e-np.cos(E))


    if deg: return a, e, i*r2d, ta*r2d, aop*r2d, raan*r2d, [year, month, day, hour]
    else: return a, e, i, ta, aop, raan, [year, month, day, hour]

def calc_epoch(epoch):

    """Fonction permettant d'extraire l'année, le mois, le jour et l'heure d'un fichier TLE

    inputs:
    epoch -- année/mois/jour/heure concaténés en une chaine de caractères  - string
    
    outputs:
    [year, month, day, hour] -- vecteur contenant l'année, le mois, le jour et l'heure sous forme de nombres entiers - 1x4
    """
    
    year = int('20'+epoch[:2])

    epoch = epoch[2:].split('.')

    day_of_year = int(epoch[0])-1
    hour = float('0.'+epoch[1])*24.0
    date = datetime.date(year, 1, 1)+datetime.timedelta(day_of_year)

    month = float(date.month)
    day = float(date.day)

    return [year, month, day, hour]

def tle2rv(tle_filename):

    """Calcul les vecteurs position et vitesse d'un fichier TLE

    inputs:
    tle_filename -- chemin relatif et nom du fichier - string
    
    outputs:
    r -- vecteur position - array 1x3 [km]
    v -- vecteur vitesse - array 1x3 [km/s]
    date -- vecteur contenant la date  - array 1x4 [jours, heures, minutes, secondes]
    """
    
    return coes2rv(tle2coes(tle_filename))
    
def plot_n_orbits(rs,
                  labels,
                  cb=pd.earth,
                  title='Figure',
                  show_plot=False,
                  save_plot=False,
                  return_ax=False):
                              
    """trace plusieurs orbites sur un même graphe 3D

    inputs:
    rs -- array contenant n vecteurs positions - array of arrays nx3 [km]
    labels -- array contenant les labels respectifs de chacune des orbites à tracer - array of strings 1xn 
    title -- titre de la figure - string
    show_plot -- impose le tracé à la fin de l'execution de la fonction - boolean (par défaut : False)
    save_plot -- impose l'enregistrement de la figure à la fin de l'execution de la fonction - boolean (par défaut : False)
    return_ax -- impose lerenvoie de l'objet "axe" à la fin de l'execution de la fonction - boolean (par défaut : False)
    
    outputs:
    (facultatif) ax -- objet matplotlib Axe
    """

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
    if return_ax:
        return ax

if __name__ == "__main__":
    r0 = np.array([-2384.46,5729.01,3050.46])
    v0 = np.array([-7.36138,-2.98997,1.64354])
    state0=np.array(rv2coes(r0,v0,print_results=True)+[[0,0,0]])