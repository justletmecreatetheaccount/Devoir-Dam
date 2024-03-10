import math
from types import prepare_class
from matplotlib.cbook import print_cycles
from scipy.integrate import odeint as odeint
import numpy as np
import matplotlib.pyplot as plt

"""
UNIQUEMENT DES GRAMMES
"""


"""
MESURES SUR LE MOTEUR
"""
RAYON_VILBREQUIN = 0.045 #m
LONGUEUR_BIELLE = 0.18 #m
DIAMETRE_CYLINDRE = 0.08 #m
VOLUME_MINIMUN = 0.0000452389 #m³

MASSE_PISTON = -1
MASSE_BIELLE = -1
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
RAPPORT_AIR_ESSENCE = 15.4 #g_air/g_essence (14.5 normalement)
MASSE_MOLAIRE_DIESEL = 142.3 #g/mole
RAPPORT_AIR_DIESEL = 26 #g_air/g_essence

VOLUME_MAX = RAYON_VILBREQUIN * 2 * math.pi * math.pow(DIAMETRE_CYLINDRE / 2, 2) + VOLUME_MINIMUN
VOLUME_COURSE = VOLUME_MAX - VOLUME_MINIMUN
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

    if (theta > -2*np.pi) and (theta < -np.pi):  #Admission
        dPression = 0
        return(dPression)
    elif (theta < thetaC):                        #Compression
        dPression = -GAMMA*pression/volume * deltaVolume
        return(dPression)
    elif (theta < thetaC + deltaThetaC):      #Combustion
        dPression = -GAMMA*pression/volume * deltaVolume + (GAMMA-1)/volume * ddQ
        return(dPression)
    elif (theta < np.pi):                        #Détente
        dPression = -GAMMA*pression/volume * deltaVolume
        return(dPression)
    else:                      #Échappement
        dPression = 10**5 - pression
        return(dPression)


        

def Pression(theta, thetaC, deltaThetaC, rapport_air_carburant):
    return odeint(Evolution_pression, PRESSION_ADMISSION, theta, args = (thetaC, deltaThetaC, rapport_air_carburant))


def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    #VOTRE CODE
    vitesse_angulaire = -1 #omega -> d theta/ d temps 
    V_output = Volume_Instant(theta)
    Q_output = Apport_Chaleur_Instant(theta, thetaC, deltaThetaC, RAPPORT_AIR_ESSENCE)
    F_pied_output = -1
    F_tete_output = -1
    p_output = np.ravel(Pression(theta, thetaC, deltaThetaC, RAPPORT_AIR_ESSENCE))
    t = -1
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


theta = np.linspace(-2*np.pi, 2*np.pi, 1001)
volume, chaleur, F, Fo, pression, t = myfunc(2555, 1.9, theta, np.radians(-26), np.radians(43))

plot_v(theta, volume)
plot_c(theta, chaleur)
plot_p(theta, pression)
