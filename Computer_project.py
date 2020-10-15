"""
Preperations
"""
import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
from sympy import *
import math
#The following two path extensions are located in the same directory as the main-file.
sys.path.append('radial.py')
sys.path.append('radiallog.py')
import radial
import radiallog


"""
Task One
"""
Z = 1
r = np.linspace(0, 35, 35 * 120)
P1s = 2 * r * np.exp(-r * Z) * Z ** (3/2)
P2s = 1/math.sqrt(2) * r * np.exp(-r * Z/2) * (1 - 1/2 * Z*r)
P2p = 1/(2 * math.sqrt(6)) * r ** 2 * np.exp(-r * Z/2) * Z ** (5/2)
P3s = 2/(3 * math.sqrt(3)) * Z ** (3/2) * r * np.exp(-r * Z / 3) * (1 - 2/3 * Z * r + 2/27 * Z ** 2 * r ** 2)
P3p = 8/(27 * math.sqrt(6)) * Z ** (5/3) * r ** 2 * np.exp(-r * Z/3) * (1 - 1/6 * Z * r)
P3d = 4/(81 * math.sqrt(30)) * Z ** (7/2) * r ** 3  * np.exp(-r * Z/3)
plt.plot(r, P1s, 'b', label = 'P1s state')
plt.plot(r, P2s, 'b', label = 'P2s state')
plt.plot(r, P2p, 'r', label = 'P2p state')
plt.plot(r, P3s, 'b', label = 'P3s state')
plt.plot(r, P3p, 'r', label = 'P3p state')
plt.plot(r, P3d, 'g', label = 'P3d State')
plt.title('Distribution for given state')
plt.ylabel(r'$P(r)$')
plt.xlabel(r'$r$')
plt.legend()
plt.show()

"""
Task 2
"""
z = (1,2,3,4)
r = np.linspace(0, 10, 10 * 120)
for Z in z:
    P1s =  2 * r * np.exp(-r * Z) * Z ** (3/2)
    plt.plot(r, P1s, label = f'P1s state for {Z}')
plt.legend()
plt.title(r'Distrubution for Z = (1,2,3,4).')
plt.ylabel(r'$P(r)$')
plt.xlabel(r'$r$')
plt.show()

"""
Task 3
"""
r = np.linspace(0, 35, 35 * 120)
Z = 1
RadialDensity_P3s = abs(2/(3 * math.sqrt(3)) * Z ** (3/2) * r * np.exp(-r * Z / 3) * (1 - 2/3 * Z * r + 2/27 * Z ** 2 * r ** 2)) ** 2
RadialDensity_P3p = abs(8/(27 * math.sqrt(6)) * Z ** (5/3) * r ** 2 * np.exp(-r * Z/3) * (1 - 1/6 * Z * r)) ** 2
RadialDensity_P3d = abs(4/(81 * math.sqrt(30)) * Z ** (7/2) * r ** 3  * np.exp(-r * Z/3)) ** 2

plt.plot(r, RadialDensity_P3s, 'b', label = 'P3s state')
plt.plot(r, RadialDensity_P3p, 'r', label = 'P3p state')
plt.plot(r, RadialDensity_P3d, 'g', label = 'P3d state')
plt.title('Probability distribution')
plt.xlabel(r'$r$')
plt.ylabel(r'$D(r) = |P(r)|^2$')
plt.legend()
plt.show()


"""
Task 4
"""
#bash -> pip install sympy

"""
Task 5
"""


R, theta, phi = symbols('R theta phi')
init_printing(use_unicode = True)
SphericalHarmonics = {
    'Y00' : 1/math.sqrt(4 * math.pi),
    'Y11' : math.sqrt(3/(8 * math.pi)) * sin(theta) * exp(-I * phi),
    'Y10' : math.sqrt(3/(4 * math.pi)) * cos(theta),
    'Y22' : math.sqrt(15/(32 * math.pi)) * sin(theta) ** 2 * exp(-2 * I * phi),
    'Y21' : math.sqrt(15/(8 * math.pi)) * sin(theta) * cos(theta) * exp(-I * phi),
    'Y20' : math.sqrt(5/(16 * math.pi)) * (2*cos(theta) ** 2 - sin(theta) ** 2)
}

Yoo_Y11 = integrate((1/(math.sqrt(4 * math.pi)) * math.sqrt(3/(8 * math.pi))*sin(theta)*exp(-I * phi)) * sin(theta), (phi, 0, 2 * math.pi), (theta, 0 , math.pi))
print(re(Yoo_Y11))
Yoo_Yoo = integrate(1/math.sqrt(4 * math.pi) * 1/math.sqrt(4 * math.pi) * sin(theta), (phi, 0, 2 * math.pi), (theta, 0 , math.pi))
print(Yoo_Yoo)

"""
for keys in SphericalHarmonics:
    print(SphericalHarmonics[keys])
"""
