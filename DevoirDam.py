import math
from types import prepare_class
from matplotlib.cbook import print_cycles
from scipy.integrate import odeint
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

"""
UNIQUEMENT DES GRAMMES
"""
D = 0.08 #[m]
C = 0.045 * 2 #[m]
L = 0.18 #[m]
mpiston = 0.8 #[kg]
mbielle = 0.9 #[kg]
tau = 11 #[-]
Mair_carb = 14.5 #[kg_air/kg_fuel]

"""
MESURES SUR LE MOTEUR
"""

RAYON_VILBREQUIN = C/2 #m
LONGUEUR_BIELLE = L #m
DIAMETRE_CYLINDRE = D #m
VOLUME_COURSE = RAYON_VILBREQUIN * 2 * math.pi * math.pow(DIAMETRE_CYLINDRE / 2, 2)

VOLUME_MINIMUN =  VOLUME_COURSE / (tau - 1)#m³

MASSE_PISTON = mpiston #Kg
MASSE_BIELLE = mbielle #Kg


"""
MESURES SUR LE MOTEUR
"""
"""
CONSTANTES
"""
SURALIMENTATION = 1.9
PRESSION_ADMISSION = SURALIMENTATION * 100000 #pascal
MASSE_MOLAIRE_AIR = 29 #g/mol
CHALEUR_MASSIQUE_AIR = 0.05414 #1/(g * K)

TEMPERATURE_EXT = 273.15 + 30 #K
CONSTANTE_GAZ_PFTS = 8.314 #1/(mol * K)
GAMMA = 1.3

POUVOIR_CALORIFIQUE = 43000 #J/g_essence (43000 normalement)
MASSE_MOLAIRE_ESSENCE = 114 #g/mole
RAPPORT_AIR_ESSENCE = 14.5 #g_air/g_essence (14.5 normalement)
MASSE_MOLAIRE_DIESEL = 142.3 #g/mole
RAPPORT_AIR_DIESEL = 26 #g_air/g_essence

VOLUME_MAX = VOLUME_COURSE + VOLUME_MINIMUN

"""
CONSTANTES
"""

def Apport_Chaleur_Instant(theta, thetaC, deltaThetaC, rapport_air_carburant):

    moles_air = (PRESSION_ADMISSION * (VOLUME_MAX - VOLUME_MINIMUN)) / (TEMPERATURE_EXT * CONSTANTE_GAZ_PFTS)
    masse_air = moles_air * MASSE_MOLAIRE_AIR #GRAMMES
    masse_carburant = masse_air / rapport_air_carburant #GRAMMES
    Q_tot = masse_carburant * POUVOIR_CALORIFIQUE
    return Q_tot * (1/2) * (1 - np.cos(math.pi * ((theta - thetaC) / deltaThetaC)))

def Evolution_apport_chaleur(theta, thetaC, deltaThetaC, rapport_air_carburant):
    moles_air = (PRESSION_ADMISSION * (VOLUME_MAX - VOLUME_MINIMUN)) / (TEMPERATURE_EXT * CONSTANTE_GAZ_PFTS)
    masse_air = moles_air * MASSE_MOLAIRE_AIR #GRAMMES
    masse_carburant = masse_air / rapport_air_carburant #GRAMMES
    Q_tot = masse_carburant * POUVOIR_CALORIFIQUE
    return (Q_tot* (np.pi*np.sin(( (theta-thetaC)*np.pi) /deltaThetaC))) / (2*deltaThetaC)


def Volume_Instant(theta):

    hauteur_piston = np.cos(theta) * RAYON_VILBREQUIN + np.cos(np.arcsin((np.sin(theta)*RAYON_VILBREQUIN)/LONGUEUR_BIELLE)) * LONGUEUR_BIELLE
    hauteur_piston_max = LONGUEUR_BIELLE + RAYON_VILBREQUIN
    dh = hauteur_piston_max - hauteur_piston
    return VOLUME_MINIMUN + math.pi * (DIAMETRE_CYLINDRE / 2)**2 * dh

def Evolution_pression(pression, theta, thetaC, deltaThetaC, rapport_air_carburant):
    ddQ = Evolution_apport_chaleur(theta, thetaC, deltaThetaC, rapport_air_carburant)
    volume = Volume_Instant(theta)
    deltaVolume = (VOLUME_COURSE*(np.sin(theta)+(np.sin(theta)*np.cos(theta))/((LONGUEUR_BIELLE/RAYON_VILBREQUIN)**2-np.sin(theta)*np.sin(theta))**0.5))/2
    dPression = 0 

    if (theta < thetaC):                        #Compression
        dPression = -GAMMA*pression/volume * deltaVolume
        return(dPression)
    elif (theta < thetaC + deltaThetaC):      #Combustion
        dPression = -GAMMA*pression/volume * deltaVolume + (GAMMA-1)/volume * ddQ
        return(dPression)
    else:                        #Détente
        dPression = -GAMMA*pression/volume * deltaVolume
        return(dPression)


def Pression(theta, thetaC, deltaThetaC, rapport_air_carburant):
    return odeint(Evolution_pression, PRESSION_ADMISSION, theta, args = (thetaC, deltaThetaC, rapport_air_carburant))

def Rankine(t, f_crit, constantes_moment, k):
    F_euler = (math.pi**2 * 200 * 10**9 * constantes_moment * t**4)/((k * LONGUEUR_BIELLE)**2)
    return 1/(F_euler) + 1/(11 * t**2 * 450 * 10**6) - 1/(f_crit)

def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    global SURALIMENTATION
    global PRESSION_ADMISSION
    SURALIMENTATION = s
    PRESSION_ADMISSION = SURALIMENTATION * 100000
    theta = np.radians(theta)
    thetaC = np.radians(-thetaC)
    deltaThetaC = np.radians(deltaThetaC)
    #VOTRE CODE
    vitesse_angulaire = rpm / 60 * 2 * math.pi #omega -> d theta/ d temps 
    V_output = Volume_Instant(theta)
    Q_output = Apport_Chaleur_Instant(theta, thetaC, deltaThetaC, Mair_carb)
    for i in range(len(theta)):
        if (theta[i] < thetaC) or (theta[i] > thetaC + deltaThetaC):
            Q_output[i] = 0
    p_output = np.ravel(Pression(theta, thetaC, deltaThetaC, Mair_carb))
    F_pression = p_output * np.pi * DIAMETRE_CYLINDRE ** 2 / 4
    F_pied_output = F_pression - MASSE_PISTON * RAYON_VILBREQUIN * vitesse_angulaire ** 2 * np.cos(theta)
    F_tete_output = -F_pression + (MASSE_PISTON + MASSE_BIELLE) * RAYON_VILBREQUIN * vitesse_angulaire ** 2 * np.cos(theta)

    F_crit = max(np.max(np.abs(F_tete_output)), np.max(np.abs(F_pied_output)))
    txx = fsolve(Rankine, 0.01,args=(F_crit, 419/12, 0.05))
    tyy = fsolve(Rankine, 0.01,args=(F_crit, 131/12, 1))
    t = min(abs(txx[0]), abs(tyy[0]))
    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t) ;

def plot_p(theta, P):
    #print(P)
    fig, axs = plt.subplots()
    plt.title("pression dans le cylindre en fonction de l'angle de vilebrequin")
    axs.plot(theta, P * (10**-5),'limegreen')
    axs.set_xlabel('Angle [radian]')
    axs.set_ylabel('Pression [bar]')
    axs.grid(True)
   # plt.savefig('Pression' + '.png')
    plt.show()


def plot_v(theta, v):
    #print(P)
    fig, axs = plt.subplots()
    plt.title("volume \"du cylindre\" en fonction de l'angle de vilebrequin")
    axs.plot(theta, v,'red')
    axs.set_xlabel('Angle [radian]')
    axs.set_ylabel('Volulme [m³]')
    axs.grid(True)
   # plt.savefig('Pression' + '.png')
    plt.show()


def plot_c(theta, c):
    #print(P)
    fig, axs = plt.subplots()
    plt.title("Q")
    axs.plot(theta, c,'blue')
    axs.set_xlabel('Angle [radian]')
    axs.set_ylabel('chaleur')
    axs.grid(True)
   # plt.savefig('Pression' + '.png')
    plt.show()

def plot_f(theta, f, f2):
    #print(P)
    fig, axs = plt.subplots()
    plt.title("f")
    axs.plot(theta, f2,'blue')
    axs.plot(theta, f,'red')
    axs.set_xlabel('Angle [radian]')
    axs.set_ylabel('newton')
    axs.grid(True)
   # plt.savefig('Pression' + '.png')
    plt.show()

theta = np.linspace(-180, 180, 1001)
volume, chaleur, ft, fp, pression, t = myfunc(2522, 1.4861300018243662, theta, 36, 45)

plot_v(theta, volume)
plot_c(theta, chaleur)
plot_p(theta, pression)
plot_f(theta, ft, fp)
print(t)
