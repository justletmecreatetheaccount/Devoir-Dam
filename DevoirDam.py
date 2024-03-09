import math
import numpy as np
import matplotlib.pyplot as plt

"""
UNIQUEMENT DES GRAMMES
"""

theta = np.linspace(-2*np.pi, 2*np.pi, 1001)

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

POUVOIR_CALORIFIQUE = 2800 #J/g_essence
MASSE_MOLAIRE_ESSENCE = 114 #g/mole
RAPPORT_AIR_ESSENCE = 14.5 #g_air/g_essence
MASSE_MOLAIRE_DIESEL = 142.3 #g/mole
RAPPORT_AIR_DIESEL = 26 #g_air/g_essence

VOLUME_MAX = RAYON_VILBREQUIN * 2 * math.pi * math.pow(DIAMETRE_CYLINDRE / 2, 2) + VOLUME_MINIMUN
"""
CONSTANTES
"""



def Apport_Chaleur_Instant(theta, thetaC, deltaThetaC, pouvoir_calorifique, rapport_air_carburant,
                           pression_admission, temperature_admission, volume_max):

    moles_air = (pression_admission * volume_max) / (temperature_admission * CONSTANTE_GAZ_PFTS)
    masse_air = moles_air * MASSE_MOLAIRE_AIR #GRAMMES
    masse_carburant = masse_air / rapport_air_carburant #GRAMMES
    Q_tot = masse_carburant * pouvoir_calorifique
    return Q_tot * (1/2) * (1 - math.cos(math.pi * ((theta - thetaC) / deltaThetaC)))



def Volume_Instant(theta, longueur_bielle, rayon_vilbrequin, diametre_cylindre, volume_minimum):

    hauteur_piston = math.cos(theta) * rayon_vilbrequin + math.cos(math.asin(rayon_vilbrequin/longueur_bielle)) * longueur_bielle
    hauteur_piston_max = longueur_bielle + rayon_vilbrequin
    dh = hauteur_piston_max - hauteur_piston
    return volume_minimum + math.pi * math.pow(diametre_cylindre / 2, 2) * dh



def Pression_Instant(volume, chaleur, chaleur_massique_air, rapport_air_carburant,
                     masse_molaire_carburant, pression_admission, temperature_admission, volume_max):

    moles_air = (pression_admission * volume_max) / (temperature_admission * CONSTANTE_GAZ_PFTS)
    masse_air = moles_air * MASSE_MOLAIRE_AIR #GRAMMES
    masse_carburant = masse_air / rapport_air_carburant #GRAMMES
    moles_carburant = masse_carburant / masse_molaire_carburant
    temperature = temperature_admission + ( chaleur / ( (masse_carburant + masse_air) * chaleur_massique_air) )
    return ((moles_air) * CONSTANTE_GAZ_PFTS * temperature) / volume
    #return ((moles_air + moles_carburant) * CONSTANTE_GAZ_PFTS * temperature) / volume



def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    #VOTRE CODE
    vitesse_angulaire = -1 #omega -> d theta/ d temps 
    V_output = Volume_Instant(theta, LONGUEUR_BIELLE, RAYON_VILBREQUIN, DIAMETRE_CYLINDRE, VOLUME_MINIMUN)
    Q_output = Apport_Chaleur_Instant(theta, thetaC, deltaThetaC, POUVOIR_CALORIFIQUE, RAPPORT_AIR_ESSENCE,
                                      TEMPERATURE_EXT, PRESSION_ADMISSION, VOLUME_MAX)
    F_pied_output = -1
    F_tete_output = -1
    p_output = Pression_Instant(V_output, Q_output, CHALEUR_MASSIQUE_AIR, RAPPORT_AIR_ESSENCE, MASSE_MOLAIRE_ESSENCE,
                                PRESSION_ADMISSION, TEMPERATURE_EXT, VOLUME_MAX)
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

pression = np.zeros_like(theta)
volume = np.zeros_like(theta)
chaleur = np.zeros_like(theta)

for i in range(len(theta)):
    V, ch, F, F, pr, t = myfunc(2555, 1.9, theta[i], 26, 43)
    volume[i] = V
    pression[i] = pr
    chaleur[i] = ch

plot_p(theta, pression)
plot_v(theta, volume)
