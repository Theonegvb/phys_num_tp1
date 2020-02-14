import matplotlib.pyplot as plt
import numpy as np


def protonSpeed(protonEnergy):
    mProton = 938.27208816
    protonSpeed = 2.998e+10 * np.sqrt(1 - (protonEnergy / mProton + 1) ** -2)
    return protonSpeed


def beta(protonSpeed):
    return protonSpeed / 2.998e+10


def gamma(protonEnergy):
    mProton = 938.27208816
    return (protonEnergy / mProton) + 1


def tMax(protonSpeed):
    mElectron = 0.51099895000
    mProton = 938.27208816
    g = gamma(protonSpeed)
    return (2 * mElectron * (g - 1)) / (
                1 + 2 * g * (mElectron / mProton) + (mElectron / mProton) ** 2)


def stoppingPower(nElectron, protonEnergy, i):
    rElectron = 2.81794033e-15 * 100
    mElectron = 0.51099895000
    pSpeed = protonSpeed(protonEnergy)
    g = gamma(pSpeed)
    b = beta(pSpeed)
    return 2 * np.pi * rElectron ** 2 * mElectron * nElectron * (1 / (b ** 2)) * \
           (np.log((2 * mElectron * (b ** 2) * (g ** 2) * tMax(pSpeed)) / (i ** 2)) - 2 * (b ** 2))


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