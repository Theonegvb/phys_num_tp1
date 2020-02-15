import scipy as scp
import numpy as np
from numpy import pi, random
from scipy import integrate
import matplotlib.pyplot as plt


def tMax(g):
    mElectron = 0.51099895000
    mPronton = 938.272
    return (2 * mElectron * ((g ** 2) - 1)) / (1 + 2 * g * (mElectron / mPronton) + (mElectron / mPronton) ** 2)


def beta(g):
    return np.sqrt(1 - 1 / (g ** 2))


def gamma(t):
    mProton = 938.272
    return t / mProton + 1


def f(x):
    i = 75 / 1000000
    rElectron = 2.81e-13
    nElectron = 3.3428e23
    mElectron = 0.51099895000

    g = gamma(x)
    b = beta(g)
    return 1 / (2 * pi * rElectron ** 2 * mElectron * nElectron * b ** -2 *
                (np.log((2 * mElectron * b ** 2 * g ** 2 * tMax(g))/(i ** 2)) - 2 * b ** 2))


energies = list(np.linspace(75, 0, 10000))

depot = np.array([])
distanceTotale = np.array([])

for i in range(len(energies) - 1):
    pas = integrate.quad(f, energies[i + 1], energies[i])[0]
    depot = np.append(depot, (energies[i] - energies[i + 1]) / pas)
    distanceTotale = np.append(distanceTotale, pas)
distanceTotale = list(np.cumsum(distanceTotale))

plt.plot(distanceTotale, depot)
plt.title("Dépot d'énergie en fonction de la profondeur pour 75 MeV.")
plt.axis(xmin=0, ymin=0, ymax=1000)
plt.show()