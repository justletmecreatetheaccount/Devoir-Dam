import math
#import matplotlib.pyplot as pyplot

"""
UNIQUEMENT DES GRAMMES
"""

"""
MESURES SUR LE MOTEUR
"""
RAYON_VILBREQUIN = -1 #m
LONGUEUR_BIELLE = -1 #m
DIAMETRE_CYLINDRE = -1 #m
VOLUME_MINIMUN = -1 #m³

MASSE_PISTON = -1
MASSE_BIELLE = -1
"""
MESURES SUR LE MOTEUR
"""
"""
CONSTANTES
"""
SURALIMENTATION = -1
PRESSION_ADMISSION = SURALIMENTATION * 100000 #pascal
MASSE_MOLAIRE_AIR = 29 #g/mol
CHALEUR_MASSIQUE_AIR = 0.05414 #1/(g * K)

TEMPERATURE_EXT = 273.15 + 25 #K
CONSTANTE_GAZ_PFTS = 8.314 #1/(mol * K)

POUVOIR_CALORIFIQUE = 43000 #J/g_essence
MASSE_MOLAIRE_ESSENCE = 114 #g/mole
RAPPORT_AIR_ESSENCE = 14.5 #g_air/g_essence
MASSE_MOLAIRE_DIESEL = 142.3 #g/mole
RAPPORT_AIR_DIESEL = 26 #g_air/g_essence

VOLUME_MAX = RAYON_VILBREQUIN * 2 * math.pi * math.pow(2, DIAMETRE_CYLINDRE / 2) + VOLUME_MINIMUN
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

    hauteur_cylindre = math.cos(theta) * (longueur_bielle + rayon_vilbrequin)
    hauteur_cylindre_max = longueur_bielle + rayon_vilbrequin
    dh = hauteur_cylindre_max - hauteur_cylindre
    return volume_minimum + math.pi * pow(2, diametre_cylindre / 2) * dh



def Pression_Instant(volume, chaleur, chaleur_massique_air, rapport_air_carburant,
                     masse_molaire_carburant, pression_admission, temperature_admission, volume_max):

    moles_air = (pression_admission * volume_max) / (temperature_admission * CONSTANTE_GAZ_PFTS)
    masse_air = moles_air * MASSE_MOLAIRE_AIR #GRAMMES
    masse_carburant = masse_air / rapport_air_carburant #GRAMMES
    moles_carburant = masse_carburant / masse_molaire_carburant
    temperature = temperature_admission + ( chaleur / ( (masse_carburant + masse_air) * chaleur_massique_air) )
    return ((moles_air + moles_carburant) * CONSTANTE_GAZ_PFTS * temperature) / volume



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
