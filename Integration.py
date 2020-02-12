import numpy as np

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
            self.array_integration = np.linspace(self.bornes[0], self.bornes[1], self.nombre_de_pas)
        else:
            self.array_integration = intervalle_variable

    def update_array_integration(self):
        """
        Fonction qui permet de mettre à jour le array qui permet de faire l'intégration.
        """
        self.array_integration = np.linspace(self.bornes[0], self.bornes[1], self.nombre_de_pas)

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
        self.nombre_de_pas = nouveau_nombre_de_pas
        self.update_array_integration()

    def get_intervalle_variable(self):
        return self.intervalle_variable

    def set_intervale_variable(self, nouvel_intervalle_variable):
        self.intervalle_variable = nouvel_intervalle_variable
        self.update_array_integration()

    def aire_rectangle(self, x_gauche, x_droite):
        return (x_droite - x_gauche) * (self.fonction_integrer(x_droite))

    def aire_triangle(self, x_gauche, x_droite):
        return (x_droite - x_gauche) * (self.fonction_integrer(x_gauche) - self.fonction_integrer(x_droite)) / 2

    def aire_trapeze (self, x_gauche, x_droite):
        return self.aire_rectangle(x_gauche, x_droite) + self.aire_triangle(x_gauche, x_droite)

    def aire_totale(self):
        aire_totale = 0
        for i in range(len(self.array_integration) - 1):
            aire_totale += self.aire_trapeze(self.array_integration[i], self.array_integration[i+1])
        return aire_totale

    def incertitude_pratique(self):
        aire_N = self.aire_totale()
        self.set_nombre_de_pas(int(2 * self.nombre_de_pas))
        aire_2N = self.aire_totale()
        self.set_nombre_de_pas(int(self.nombre_de_pas / 2))

        return (aire_2N - aire_N)/3




def fonction_test(x):
    return x ** 2 - x ** 3
i=5
while True:
    integration_1 = Integration_trapezes(fonction_test, bornes=[0, np.pi], nombre_de_pas=i)
    if integration_1.incertitude_pratique() <= 0.01:
        break
    i += 1

print(integration_1.nombre_de_pas)
print(integration_1.aire_totale())
