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
ss = np.array([])
for i in range(len(energies) - 1):
    s = integrate.quad(f, energies[i], energies[i + 1])[0]
    ss = np.append(ss, s)
ss = list(np.cumsum(ss))
energies.pop(0)

plt.plot(energies, ss)
plt.show()