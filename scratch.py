import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

#%matplotlib inline


def beta(g):
    return np.sqrt(1 - 1 / (g ** 2))


def gamma(t):
    mProton = 938.272
    return t / mProton + 1


def tMax(g):
    mElectron = 0.51099895000
    mPronton = 938.272
    return (2 * mElectron * ((g ** 2) - 1)) / (1 + 2 * g * (mElectron / mPronton) + (mElectron / mPronton) ** 2)


def stoppingPower(nElectron, protonEnergy, i):
    rElectron = 2.81794033e-13
    mElectron = 0.51099895000
    g = gamma(protonEnergy)
    b = beta(g)
    return 2 * pi * rElectron ** 2 * mElectron * nElectron * b ** -2 *\
           (np.log((2 * mElectron * b ** 2 * g ** 2 * tMax(g))/(i ** 2)) - 2 * b ** 2)


if __name__ == '__main__':
    energies = np.arange(0.001, 10000, 0.0005)

    bHydrogen = (1 * 0.047234 * 6.022e23) / (1.660540199e-24 * 6.022e23)
    bCarbon = (6 * 0.144330 * 6.022e23) / (1.994242358e-23 * 6.022e23)
    bNitrogen = (7 * 0.041990 * 6.022e23) / (2.325824007e-23 * 6.022e23)
    bOxygen = (8 * 0.446096 * 6.022e23) / (2.656703247e-23 * 6.022e23)
    bMagnesium = (12 * 0.002200 * 6.022e23) / (4.035776902e-23 * 6.022e23)
    bPhosphorus = (15 * 0.104970 * 6.022e23) / (5.143317694e-23 * 6.022e23)
    bSulfur = (16 * 0.003150 * 6.022e23) / (5.323525827e-23 * 6.022e23)
    bCalcium = (20 * 0.209930 * 6.022e23) / (6.655113013e-23 * 6.022e23)
    bZinc = (30 * 0.000100 * 6.022e23) / (1.085661182e-22 * 6.022e23)
    bNElectron = 1.85 * (bHydrogen + bCarbon + bNitrogen + bOxygen + bMagnesium + bPhosphorus + bSulfur + bCalcium + bZinc)
    bI = 106.4 / 1000000
    plt.plot(energies, stoppingPower(bNElectron, energies, bI))
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

    hydrogen = (1 * 0.111894 * 6.022e23) / (1.660540199e-24 * 6.022e23)
    oxygen = (8 * 0.888106 * 6.022e23) / (2.656703247e-23 * 6.022e23)
    nElectron = 1 * (hydrogen + oxygen)
    i = 75 / 1000000
    plt.plot(energies, stoppingPower(nElectron, energies, i))
    plt.xscale('log')
    plt.yscale('log')
    plt.show()