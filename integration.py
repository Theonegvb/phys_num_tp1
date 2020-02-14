# Tests
import numpy as np
import matplotlib.pyplot as plt
import scipy as scp
from numpy import linspace, cos, pi, tan, ones, copy, random

class Integration:
    """
    Classe qui permet de faire l'intégration de fonctions réelle.
    """

    def __init__(self, fonction_integrer=1, bornes=[1,2], nombre_de_pas=1, intervalle_variable=None):
        """
        Le constructeur de l'intégration.
        Args:
            fonction_integrer (fonction): La fonction qui sera à intégrer numériquement.
            bornes (array like): Donne les bornes d'intégration numérique.
            nombre_de_pas (int): Nombre de pas dans l'intégration numérique.
            emplacement_variable (array like): On peut donné directement l'emplacement des pas dans l'inégration entre
             les bornes.
        """

        self.intervalle_variable = intervalle_variable
        self.fonction_integrer=fonction_integrer
        self.bornes = bornes
        self.nombre_de_pas = nombre_de_pas



class Integration_trapezes(Integration):
    """
    Classe de l'intégration par méthode des trapèzes
    """
    def __init__(self, fonction_integrer=1, bornes=np.array([1, 2]), nombre_de_pas=1, intervalle_variable=None):
        """
        Le constructeur de l'intégration.
        Args:
            fonction_integrer (fonction): La fonction qui sera à intégrer numériquement.
            bornes (array like): Donne les bornes d'intégration numérique.
            nombre_de_pas (int): Nombre de pas dans l'intégration numérique.
            emplacement_variable (array like): On peut donné directement l'emplacement des pas dans l'inégration entre
             les bornes.
        """
        super().__init__(fonction_integrer, bornes, nombre_de_pas, intervalle_variable)
        if intervalle_variable == None:
            self.array_integration = np.linspace(self.bornes[0], self.bornes[1], self.nombre_de_pas + 1)
        else:
            self.array_integration = intervalle_variable

        array_valeur_fonction = np.array([])
        for i in self.array_integration:
            array_valeur_fonction = np.append(array_valeur_fonction, self.fonction_integrer(i))

        self.array_valeur_fonction = array_valeur_fonction

    def update_array_integration(self):
        self.array_integration = np.linspace(self.bornes[0], self.bornes[1], self.nombre_de_pas + 1)

    def update_array_valeur_fonction(self):
        array_valeur_fonction = np.array([])
        for i in self.array_integration:
            array_valeur_fonction = np.append(array_valeur_fonction, self.fonction_integrer(i))

        self.array_valeur_fonction = array_valeur_fonction

    def get_array_valeur_fonction(self):
        return self.array_valeur_fonction

    def get_array_integration(self):
        return self.array_integration

    def get_fonction_integrer(self):
        return self.fonction_integrer

    def set_fonction_integrer(self, nouvelle_fonction_integrer):
        self.fonction_integrer = nouvelle_fonction_integrer
        self.update_array_integration()

    def get_bornes(self):
        return self.bornes

    def set_bornes(self, nouvelles_bornes):
        self.bornes = nouvelles_bornes
        self.update_array_integration()

    def get_nombre_de_pas(self):
        return self.nombre_de_pas

    def set_nombre_de_pas(self, nouveau_nombre_de_pas):

        if nouveau_nombre_de_pas/self.nombre_de_pas % 2 == 0:
            nombre_de_fois_par_2 = int(np.floor(np.log(nouveau_nombre_de_pas/self.nombre_de_pas)/np.log(2)))
            if nombre_de_fois_par_2 == np.log(nouveau_nombre_de_pas)/np.log(2):
                for i in range(nombre_de_fois_par_2):
                    self.double_nombre_de_pas()
            self.nombre_de_pas = nouveau_nombre_de_pas
        else:
            self.nombre_de_pas = nouveau_nombre_de_pas
            self.update_array_integration()
            self.update_array_valeur_fonction()

    def double_nombre_de_pas(self):
        self.nombre_de_pas = 2 * self.nombre_de_pas
        self.update_array_integration()

        nouvel_array_valeur_fonction = np.array([])

        for i in range(self.nombre_de_pas + 1):
            if i % 2 == 0:
                nouvel_array_valeur_fonction = \
                    np.append(nouvel_array_valeur_fonction, self.array_valeur_fonction[int(i/2)])
            else:
                nouvel_array_valeur_fonction = \
                    np.append(nouvel_array_valeur_fonction, self.fonction_integrer(self.array_integration[i]))

        self.array_valeur_fonction = nouvel_array_valeur_fonction

    def get_intervalle_variable(self):
        return self.intervalle_variable

    def set_intervale_variable(self, nouvel_intervalle_variable):
        self.intervalle_variable = nouvel_intervalle_variable
        self.update_array_integration()

    def aire_trapeze (self, index_droite):
        return (self.array_integration[index_droite] - self.array_integration[index_droite-1])\
               * (self.array_valeur_fonction[index_droite] + self.array_valeur_fonction[index_droite-1]) / 2

    def aire_totale(self):
        aire_totale = 0
        for index_droite in range(len(self.array_integration) - 1):
            aire_totale += self.aire_trapeze(index_droite+1)
        return aire_totale

    def incertitude_pratique(self):
        aire_N = self.aire_totale()
        self.set_nombre_de_pas(int(2 * self.nombre_de_pas))
        aire_2N = self.aire_totale()
        self.set_nombre_de_pas(int(self.nombre_de_pas / 2))

        return (aire_2N - aire_N)/3






def tMax(g):
    mElectron = 0.51099895000
    mPronton = 938.272
    return (2 * mElectron * ((g ** 2) - 1)) / (1 + 2 * g * (mElectron / mPronton) + (mElectron / mPronton) ** 2)

def beta(g):
    return np.sqrt(1 - 1 / (g ** 2))

def gamma(t):
    mProton = 938.272
    return t / mProton + 1

def S_col_eau(x):
    i = 75 / 1000000
    rElectron = 2.81e-13
    nElectron = 3.3428e23
    mElectron = 0.51099895000

    g = gamma(x)
    b = beta(g)


    return 1 / (2 * pi * rElectron ** 2 * mElectron * nElectron * b ** -2 *
                (np.log((2 * mElectron * b ** 2 * g ** 2 * tMax(g))/(i ** 2)) - 2 * b ** 2))


integration_2 = Integration_trapezes(S_col_eau, bornes=[0.001, 100], nombre_de_pas=1000)

print(integration_2.aire_totale())
