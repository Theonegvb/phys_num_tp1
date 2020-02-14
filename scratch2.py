import numpy as np
import scipy as scp
import matplotlib.pyplot as plt
from numpy import linspace, cos, pi, tan, ones, copy, random

def gaussxw(N):
    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3, 4 * N - 1, N) / (4 * N + 2)
    x = cos(pi * a + 1 / (8 * N * N * tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N, float)
        p1 = copy(x)
        for k in range(1,N):
            p0, p1 = p1,((2 * k + 1) * x * p1 - k * p0)/(k + 1)
        dp = (N + 1)*(p0 - x * p1)/(1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2 * (N + 1) * (N + 1) / (N * N * (1 - x * x) * dp * dp)

    return x, w

def gaussianProtonsGenerator():
    return random.normal(240, 3, 10000)

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

def gaussianQuadrature(N, a, b):
    x, w = gaussxw(N)
    xp = 0.5 * (b - a) * x + 0.5 * (b + a)
    wp = 0.5 * (b - a) * w

    # Perform the integration
    s = 0
    for k in range(N):
        s += wp[k] * f(xp[k])
    return s

def optimalQuadrature(a, b):
    N = 10
    s1 = gaussianQuadrature(N, a, b)
    N = 2 * N
    while True:
        s2 = gaussianQuadrature(N, a, b)
        delta = s1 - s2
        print(N, delta, s1, s2)
        if abs(delta) <= 1e-9:
            break
        N = 2 * N
        if N > 1000:
            break
    return N

if __name__ == '__main__':
    N = optimalQuadrature(0, 240)

    protons = gaussianProtonsGenerator()
    distances = []
    for proton in protons:
        N = 10
        a = 0
        b = proton
        distance = gaussianQuadrature(N, a, b)
        distances.append(distance)

    plt.hist(distances, 100)
    plt.xlabel(r'Distance [$cm$]')
    plt.legend('show')
    plt.show()