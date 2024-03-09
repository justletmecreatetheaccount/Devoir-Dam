from scipy.integrate import odeint as odeint
import numpy as np
from matplotlib import pyplot as plt

"""valeurs initiales"""

tau = 11 #@valeur taux compression@ #[-]
D = 0.08 #@valeur alesage@ #[m]
C = 0.09 #@valeur course@ #[m]
L = 0.18 #@valeur longueur bielle@ #[m]
mpiston = 0.8 #@valeur masse piston@ #[kg]
mbielle = 0.9 #@valeur masse bielle@ #[kg]
Q = 2800000 #@valeur chaleur emise par fuel par kg de melance admis@ #[J/kg_inlet gas]


"""code"""

D = D
L = L             # longueur de la bielle, en [m]
R = C/2           # longueur de la manivelle, en [m]
Vc = np.pi * D**2 * R / 2     #difference entre volume max et vol min
#Qtot = Q*0.75*Vc      # en J (chaleur émise par fuel par kg * masse volumique de l'essence * Vc /1000)
mPiston = mpiston #kg
mBielle = mbielle #kg
gamma = 1.3
beta = L/R
r = tau  

theta = np.linspace(-2*np.pi, 2*np.pi, 1001)

def Qtot(s):
    Nmoles = (s*Vc)/(8.314*303)
    Mgaz = Nmoles*0.029
    Qto = Mgaz*Q
    return Qto


def Volume(theta):
    return (Vc*(1-np.cos(theta)+beta-(beta*beta-np.sin(theta)*np.sin(theta))**0.5))/2+Vc/(r-1)


def dQ(theta, thetaC, deltaThetaC,s):
    return (Qtot(s)*(np.pi*np.sin(((theta-thetaC)*np.pi)/deltaThetaC)))/(2*deltaThetaC)


def dV(theta):
    return (Vc*(np.sin(theta)+(np.sin(theta)*np.cos(theta))/(beta*beta-np.sin(theta)*np.sin(theta))**0.5))/2


def dpdtheta(p, theta, thetaC, deltaThetaC, s):
    """
    retourne la dérivé de la pression en fonction du temps entre -2pi et +2pi
   
    """
    dq = dQ(theta, thetaC, deltaThetaC,s)
    V = Volume(theta)
    dv = dV(theta)
    dpdt = 0 

    if (theta > -2*np.pi) and (theta < -np.pi):  #Admission
        dpdt = 0
    elif (theta < thetaC):                        #Compression
        dpdt = -gamma*p/V*dv
    elif (theta < thetaC + deltaThetaC):      #Combustion
        dpdt = -gamma*p/V*dv + (gamma-1)/V*dq
    elif (theta < np.pi):                        #Détente
        dpdt = -gamma*p/V*dv
    elif (theta < 2*np.pi):                      #Échappement
        dpdt = s-p
        
    return(dpdt)

def p(theta, thetaC, deltaThetaC, s):
    """
    retourne la pression entre -2pi et +2pi en integrant dpdtheta
    """
    p = odeint(dpdtheta, s, theta, args = (thetaC, deltaThetaC, s))
    return p   

def forcemax(theta, tpm, thetaC, deltaThetaC, s):
    """
    retourne la valeur de la force maximale s'exerçant sur la bielle
    """
    omega = tpm*2*np.pi/(60)  #omega en rad/s
    pression = np.ravel(p(theta, thetaC, deltaThetaC, s))
    plot_p(theta, pression)
    
    Fp = np.pi*D**2*pression/4
    FmPiston = mPiston*R*omega**2*np.cos(theta)         
    FmBielle = mBielle*R*omega**2*np.cos(theta)
    
    Fpied = Fp - FmPiston                           
    Ftete = -Fp + FmBielle + FmPiston
    dessine(theta, Ftete, Fpied)
    
    Fpied_max = max(Fpied)
    Ftete_max = max(-Ftete)
    Fmax = max(Fpied_max, Ftete_max)
    
    return Fmax

def myfunc(rpm, s, thetaC, deltaThetaC):
    C_s = 1         #coefficient de sécurité visblement il faut pas le prendre en compte
    t = 0           # inconnue à trouver
    thetaC = np.radians(-thetaC)
    deltaThetaC = np.radians(deltaThetaC)
    s = s*10**5

    F_crit = forcemax(theta, rpm, thetaC, deltaThetaC, s)*C_s
    
    Feulerxx = (np.pi*np.pi*200*(10**9)*(419/12)*t**4)/(1*L)**2
    Feuleryy = (np.pi*np.pi*200*(10**9)*(131/12)*t**4)/(0.5*L)**2
    
    A = 11*t**2    #aire de la section
    
    Poly1 = [1/F_crit,0,-1/(11*450*10**6),0,-1/((np.pi*np.pi*200*(10**9)*(419/12))/(1*L)**2)] # équation à résoudre
    Poly2 = [1/F_crit,0,-1/(11*450*10**6),0,-1/((np.pi*np.pi*200*(10**9)*(131/12))/(0.5*L)**2)]
    Rac1 = np.roots(Poly1)
    Rac2 = np.roots(Poly2)
    
    racReal1=[]
    racReal2=[]
    
    for i in Rac1:
        if (i.imag==0): #Pour prendre seulement les racines réelles 
            racReal1.append(i.real)
    for j in Rac2:
        if(j.imag==0):
            racReal2.append(j.real)
            
    t=max([max(racReal1),max(racReal2)])
    
    return(t)

    
def dessine(theta, Ftete, Fpied):  
    fig, axs = plt.subplots()
    plt.title("efforts")
    p1 = axs.plot(theta, Ftete,'limegreen', label = 'Ftete')
    p2 = axs.plot(theta, Fpied, 'blue', label = 'Fpied')
    axs.set_xlabel('Angle [radian]')
    axs.set_ylabel('Forces[N]')
    axs.grid(True)
   # plt.savefig('Pression' + '.png')
    plt.legend()
    plt.show()

def plot_p(theta, P):
    #print(P)
    fig, axs = plt.subplots()
    plt.title("pression dans le cylindre en fonction de l'angle de vilebrequin")
    axs.plot(theta, P*10**-5,'limegreen')
    axs.set_xlabel('Angle [radian]')
    axs.set_ylabel('Pression [bar]')
    axs.grid(True)
   # plt.savefig('Pression' + '.png')
    plt.show()

print(myfunc(2555, 1.9, 26, 43))