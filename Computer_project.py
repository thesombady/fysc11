"""
Preperations
"""
import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
from sympy import *
import math
import functools
#The following two path extensions are located in the same directory as the main-file.
#sys.path.append('radial.py')
#sys.path.append('radiallog.py')
import radial
import radiallog
import time


"""
Task One
"""
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

"""
Task 2
"""
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
"""
Task 3
"""
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

"""
Task 4
"""
#bash -> pip install sympy

"""
Task 5
"""


R, theta, phi = symbols('R theta phi')
init_printing(use_unicode = True)
"""
SphericalHarmonics = {
    'Y00' : 1/sqrt(4 * math.pi),
    'Y1+1' : -sqrt(3/(8 * math.pi)) * sin(theta) * exp(+I * phi),
    'Y10' : sqrt(3/(4 * math.pi)) * cos(theta),
    'Y2+2' : sqrt(15/(32 * math.pi)) * sin(theta) ** 2 * exp(2 * I * phi),
    'Y2+1' : -sqrt(15/(8 * math.pi)) * sin(theta) * cos(theta) * exp(+I * phi),
    'Y20' : sqrt(5/(16 * math.pi)) * (2*cos(theta) ** 2 - sin(theta) ** 2)
}

for key1 in SphericalHarmonics.keys():
    for key2 in SphericalHarmonics.keys():
        try:
            inte1 = integrate(conjugate(SphericalHarmonics[key1]) * SphericalHarmonics[key2] * sin(theta), (phi, 0, 2 * math.pi), (theta, 0 , math.pi))
            print(inte1, key1, key2)
        except:
            print(f'Cannot compute {key1}, with {key2}')
            #Could implement quad from scipy to compute because of polynomial division error for some
"""

"""
Task 6
"""

#This has to be completed
#Derivation

"""
Task 7
"""
#Studied section 10, and also the python script.

"""
Task 8
"""
"""
Z = 1

Potential = lambda r,l: -Z/r + l * (l + 1)/(2*r ** 2)
def Potential(r,l):
    return -Z/r + l * (l + 1)/(2*r ** 2)

r = np.linspace(0.00001, 5, 10 * 50)

slist = list(map(functools.partial(Potential, l = 0), r))
plist = list(map(functools.partial(Potential, l = 1), r))
dlist = list(map(functools.partial(Potential, l = 2), r))
flist = list(map(functools.partial(Potential, l = 3), r))
plt.plot(r, slist, 'r', label = 'S orbital')
plt.plot(r, plist, 'b', label = 'P orbital')
plt.plot(r, dlist, 'g', label = 'D orbital')
plt.plot(r, flist, 'y', label = 'F orbital')
plt.title(r'$V(r)$ versus $r$')
plt.ylabel(r'$V(r)$')
plt.xlabel(r'$r$')
plt.ylim(-10,10)
plt.legend()
plt.show()
"""

"""
Task 9
"""
#Save all images
Z = 1

EnergyFormula = lambda Z,n: -Z ** 2 / (2 * n ** 2)
"""
Value1 = radial.radial(0, 1, Z, plot = True, showplot = True)
print(EnergyFormula(Z, 1)) #Correct energy, E approx -0.5
Value2 = radial.radial(0, 2, Z, plot = True, showplot = True)
print(EnergyFormula(Z, 2))#Correct energy, E approx -0.125
Value3 = radial.radial(0, 3, Z, plot = True, showplot = True)
print(EnergyFormula(Z, 3))#Correct energy, E approx -0.055
Value4 = radial.radial(0, 4, Z, plot = True, showplot = True)
print(EnergyFormula(Z, 4))#Correct energy, E approx -0.0312
Value5 = radial.radial(0, 6, Z, plot = True, showplot = True)
print(EnergyFormula(Z, 6))#Correct energy, E approx - 0.0138
Value6 = radial.radial(0, 9, Z, plot = True, showplot = True)
print(EnergyFormula(Z, 9))# Correct energy, E  approx -0.00617
#Number of nodes found in the Program
print(Value1[-2])
print(Value2[-2])
print(Value3[-2])
print(Value4[-2])
print(Value5[-2])
print(Value6[-2])
#Energy
print(Value1[-1])
print(Value2[-1])
print(Value3[-1])
print(Value4[-1])
print(Value5[-1])
print(Value6[-1])
"""



#Small difference in accuracy
"""
Task 10
"""
#Derivation

"""
Task 11
"""
#The program was studied in detail, radiallog

"""
Task 12
"""

"""
EnergyFormula = lambda Z,n: -Z ** 2 / (2 * n ** 2)
Value1 = radiallog.radiallog(0, 1, Z, N=1, plot = True) #Gridpoints = 657
Value2 = radiallog.radiallog(0, 2, Z, N=1, plot = True) #Gridpoints = 690
Value3 = radiallog.radiallog(0, 3, Z, N=1, plot = True) #Gridpoints = 710
Value4 = radiallog.radiallog(0, 4, Z, N=1, plot = True) #Gridpoints = 724
Value5 = radiallog.radiallog(0, 6, Z, N=1, plot = True) #Gridpoints = 743
Value6 = radiallog.radiallog(0, 9, Z, N=1, plot = True) #Gridpoints = 763
#Energies
print(Value1[-1])
print(Value2[-1])
print(Value3[-1])
print(Value4[-1])
print(Value5[-1])
print(Value6[-1])
"""

"""
Task 13
"""
Z = 1
nValues = [1,2,3,4,6,9]
"""
for n in nValues:
    r, p, e = radiallog.radiallog(0, n, Z, N = 1, plot = False)
    plt.plot(r, np.sqrt(r) * p, label = f'{n}s orbital')


plt.title("Radiallog for different states")
plt.xlim(0, 100)
plt.xlabel(r'$R$')
plt.ylabel('Radial orbit')
plt.legend()
plt.show()
"""


"""
Task 14
"""
"""
nValues = [2,3,4,5,7,9]
p = 1
Z = 1
ListValue = []

for n in nValues:
    Value = radiallog.radiallog(p, n, Z, N=1, plot = True)
    ListValue.append(Value)
    print(EnergyFormula(Z,n))

#Gridpoints = 690, Energy = -0.12500000036952086, E = - 0.125
#Gridpoints = 710, Energy = -0.55555556565547115, E = - 0.05555555555555555
#Gridpoints = 724, Energy = -0.31250001979224021, E = - 0.03125
#Gridpoints = 734, Energy = -0.20000003269382283, E = - 0.02
#Gridpoints = 750, Energy = -0.10204088402754037, E = - 0.01020408163265306
#Gridpoints = 763, Energy = -0.061728509786443397, E = -0.006172839506172839
#Radiallog seems to be more accurate
"""

"""
Task 15
"""
"""
The effective potential given by the function Potential, depends on on n,z,l. We take the derivitve with respect ot l and see why.
d(V_{eff})/dl = (2l + 1)/2r^2
"""

"""
Task 16
"""
Z = 1
RadialFunctions = {
    'P1s' : 2 * Z ** (3/2) * R * exp(-Z * R),
    'P2s' : 1/sqrt(2) * Z ** (3/2) * R * exp(- Z * R / 2) * (1 - 1/2 * Z * R),
    'P2p' : 1/(2*sqrt(6)) * Z ** (5/2) * R ** 2 * exp(- Z * R / 2),
    'P3s' : 2/(3*sqrt(3)) * Z ** (3/2) * R * exp(- Z * R / 3) * (1 - 2/3 * Z * R +2/27 * Z ** 2 * R ** 2),
    'P3p' : 8/(27 * sqrt(6)) * Z ** (5/2) * R ** 2 * exp(- Z * R/3)*(1 - 1/6 * Z * R),
    'P3d' : 4/(81 * sqrt(30)) * Z ** (7/2) * R ** 3 * exp(- Z * R/3),
}
"""
for key in RadialFunctions.keys():
    try:
        value = integrate(R * conjugate(RadialFunctions[key]) * RadialFunctions[key], (R, 0, oo ))
        print(value, key, '\n')
    except:
        try:
            value = integrate(R * (RadialFunctions[key]) ** 2, (R, 0, oo))
            try:
                print(value, key, '\n')
            except:
                raise KeyError("Cannot print")
        except:
            print(f"Can't compute {key}")
"""
"""
1.50000000000000 P1s #Correct
6.00000000000000 P2s #Correct
5.00000000000000 P2p #Correct
13.5000000000000 P3s #Correct
12.5000000000000 P3p #Correct
10.5000000000000 P3d #correct
"""
#Compare to according plot above.
#Compute integrals by hand.

"""
Task 17
"""
#Modified the radiallog.
# Numerov g function was modified to the following:
# g = -2*r**2*(E + (Z/r) - (Z-zeta)*r/(a ** 2 + r ** 2)) + l ** 2  + l + 1/4
"""
Task 18
"""

"""
Energy1s = radiallog.radiallog(0, 1, 11, 11, plot = False, updated = True)
Energy2s = radiallog.radiallog(0, 2, 11, 11, plot = False, updated = True)
Energy2p = radiallog.radiallog(1, 2, 11, 11, plot = False, updated = True)
Energy3s = radiallog.radiallog(0, 3, 11, 11, plot = False, updated = True)
Energy3p = radiallog.radiallog(1, 3, 11, 11, plot = False, updated = True)
Energy3d = radiallog.radiallog(2, 3, 11, 11, plot = False, updated = True)
Energy4s = radiallog.radiallog(0, 4, 11, 11, plot = False, updated = True)
Energy4p = radiallog.radiallog(1, 4, 11, 11, plot = False, updated = True)
Energy4d = radiallog.radiallog(2, 4, 11, 11, plot = False, updated = True)
Energy4f = radiallog.radiallog(3, 4, 11, 11, plot = False, updated = True)
Energy5s = radiallog.radiallog(0, 5, 11, 11, plot = False, updated = True)
plt.plot(Energy1s[0], Energy1s[1], 'r', label = '1s')
plt.plot(Energy2s[0], Energy2s[1], 'b', label = '2s')
plt.plot(Energy2p[0], Energy2p[1], 'y', label = '2p')
plt.plot(Energy3s[0], Energy3s[1], 'g', label = '3s')
plt.legend()
plt.xlim(0,5)
plt.ylim(-2,4)
plt.title('Task 18')
plt.ylabel(r'$P(r)$, $[au]$')
plt.xlabel(r'$r$, $[au]$')
plt.show()

print(Energy1s[-1])
print(Energy2s[-1])
print(Energy2p[-1])
print(Energy3s[-1])
print(Energy3p[-1])
print(Energy3d[-1])
print(Energy4s[-1])
print(Energy4p[-1])
print(Energy4d[-1])
print(Energy4f[-1])
print(Energy5s[-1])
"""

"""
Task 19
"""

cm1converter = 219474.63
IonizationLimit = - 229445.71 #cm^(-1
ValueOfA = []
a = 0.19
values = []
"""
while a < 0.25:
    Al1s = radiallog.radiallog(0, 1, 13, 11, plot = False, updated = True, a = a)
    Al2s = radiallog.radiallog(0, 2, 13, 11, plot = False, updated = True, a = a)
    Al2p = radiallog.radiallog(1, 2, 13, 11, plot = False, updated = True, a = a)
    Al3s = radiallog.radiallog(0, 3, 13, 11, plot = False, updated = True, a = a)
    Al3p = radiallog.radiallog(1, 3, 13, 11, plot = False, updated = True, a = a)
    Al3d = radiallog.radiallog(2, 3, 13, 11, plot = False, updated = True, a = a)
    Al4s = radiallog.radiallog(0, 4, 13, 11, plot = False, updated = True, a = a)
    Al4p = radiallog.radiallog(1, 4, 13, 11, plot = False, updated = True, a = a)
    Al4d = radiallog.radiallog(2, 4, 13, 11, plot = False, updated = True, a = a)
    Al4f = radiallog.radiallog(3, 4, 13, 11, plot = False, updated = True, a = a)
    Al5s = radiallog.radiallog(0, 5, 13, 11, plot = False, updated = True, a = a)
    EnergyOfGroundstate = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al3s[-1])
    EnergyOfExcited = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al3d[-1])
    value = abs(EnergyOfGroundstate - EnergyOfExcited) * cm1converter
    print(value - 115956.21, a)
    if abs(value - 115956.21) < 1:
        print(a)
        ValueOfA.append(a)
        break
    values.append((value,a))
    a += 0.0001
    print(value)
"""
"""
a = 0.20789999999999803
Al1s = radiallog.radiallog(0, 1, 13, 11, plot = False, updated = True, a = a)
Al2s = radiallog.radiallog(0, 2, 13, 11, plot = False, updated = True, a = a)
Al2p = radiallog.radiallog(1, 2, 13, 11, plot = False, updated = True, a = a)
Al3s = radiallog.radiallog(0, 3, 13, 11, plot = False, updated = True, a = a)
Al3p = radiallog.radiallog(1, 3, 13, 11, plot = False, updated = True, a = a)
Al3d = radiallog.radiallog(2, 3, 13, 11, plot = False, updated = True, a = a)
Al4s = radiallog.radiallog(0, 4, 13, 11, plot = False, updated = True, a = a)
Al4p = radiallog.radiallog(1, 4, 13, 11, plot = False, updated = True, a = a)
Al4d = radiallog.radiallog(2, 4, 13, 11, plot = False, updated = True, a = a)
Al4f = radiallog.radiallog(3, 4, 13, 11, plot = False, updated = True, a = a)
Al5s = radiallog.radiallog(0, 5, 13, 11, plot = False, updated = True, a = a)
EnergyOfGroundstate = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al3s[-1])
SecondExited = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al3p[-1]) * cm1converter - EnergyOfGroundstate*cm1converter
ThirdExcited = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al3d[-1]) * cm1converter - EnergyOfGroundstate*cm1converter
FourthExcited = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al4s[-1]) * cm1converter - EnergyOfGroundstate*cm1converter
FifthExcited = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al4p[-1]) * cm1converter - EnergyOfGroundstate*cm1converter
SixthExcited = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al4d[-1]) * cm1converter - EnergyOfGroundstate*cm1converter
SevenExcited = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al4f[-1]) * cm1converter - EnergyOfGroundstate*cm1converter
EightExcited = (2 * Al1s[-1] + 2 * Al2s[-1] + 6 * Al2p[-1] + 1 * Al5s[-1]) * cm1converter - EnergyOfGroundstate*cm1converter
print(SecondExited)
print(ThirdExcited)
print(FourthExcited)
print(FifthExcited)
print(SixthExcited)
print(SevenExcited)
print(EightExcited)
ListOfStates = [Al1s, Al2s, Al2p, Al3s, Al3p, Al3d, Al4s, Al4p, Al4d, Al4f, Al5s]
ListOfNames = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f', '5s']
for i in range(len(ListOfStates)):
    plt.plot(ListOfStates[i][0], ListOfStates[i][1], '-', label = ListOfNames[i])
plt.legend()
plt.title('Different States for Al')
plt.xlabel('r [a.u]')
plt.ylabel('P(r) [a.u]')
plt.xlim(0, 5)
plt.ylim(-2, 5)
plt.show()
print(Al1s[-1])
print(Al2s[-1])
print(Al2p[-1])
print(Al3s[-1])
print(Al3p[-1])
print(Al3d[-1])
print(Al4s[-1])
print(Al4p[-1])
print(Al4d[-1])
print(Al4f[-1])
print(Al5s[-1])
"""
"""
Task 20
"""
"""
a = 0.27
values = []
difference = []
while a < 0.3:
    Li1s = radiallog.radiallog(0, 1, 3, 3, plot = False, updated = True, a = a)
    Li2s = radiallog.radiallog(0, 2, 3, 3, plot = False, updated = True, a = a)
    Li2p = radiallog.radiallog(1, 2, 3, 3, plot = False, updated = True, a = a)
    Li3s = radiallog.radiallog(0, 3, 3, 3, plot = False, updated = True, a = a)
    Li3p = radiallog.radiallog(1, 3, 3, 3, plot = False, updated = True, a = a)
    Li3d = radiallog.radiallog(2, 3, 3, 3, plot = False, updated = True, a = a)
    Li4s = radiallog.radiallog(0, 4, 3, 3, plot = False, updated = True, a = a)
    Li4p = radiallog.radiallog(1, 4, 3, 3, plot = False, updated = True, a = a)
    Groundstate = 2 * Li1s[-1] + 1 * Li2s[-1]
    Excited = 2 * Li1s[-1] + 1 * Li2p[-1]
    diff = abs(Groundstate - Excited) * cm1converter
    difference.append((diff,a))
    if abs(diff - 14903.66) < 100:
        print(a)
        print('Converged')
        break
    a += 0.001
print(difference)
print(a)
Li1s = radiallog.radiallog(0, 1, 3, 3, plot = False, updated = True, a = a)
Li2s = radiallog.radiallog(0, 2, 3, 3, plot = False, updated = True, a = a)
Li2p = radiallog.radiallog(1, 2, 3, 3, plot = False, updated = True, a = a)
Li3s = radiallog.radiallog(0, 3, 3, 3, plot = False, updated = True, a = a)
Li3d = radiallog.radiallog(2, 3, 3, 3, plot = False, updated = True, a = a)
Li4s = radiallog.radiallog(0, 4, 3, 3, plot = False, updated = True, a = a)
Li4p = radiallog.radiallog(1, 4, 3, 3, plot = False, updated = True, a = a)
Groundstate = 2 * Li1s[-1] + 1 * Li2s[-1]
FirstExcited = 2 * Li1s[-1] + 1 * Li2p[-1]
SecondExcited = 2 * Li1s[-1] + 1 * Li3s[-1]
ThirdExcited = 2 * Li1s[-1] + 1 * Li3p[-1]
FourthExcited = 2 * Li1s[-1] + 1 * Li3d[-1]
FifthExcited = 2 * Li1s[-1] + 1 * Li4s[-1]
SixthExcited = 2 * Li1s[-1] + 1 * Li4p[-1]
print(FirstExcited * cm1converter - Groundstate * cm1converter)
print(SecondExcited * cm1converter - Groundstate * cm1converter)
print(ThirdExcited * cm1converter - Groundstate * cm1converter)
print(FourthExcited * cm1converter - Groundstate * cm1converter)
print(FifthExcited * cm1converter - Groundstate * cm1converter)
print(SixthExcited * cm1converter - Groundstate * cm1converter)
ListOfStates = [ Li1s, Li2s, Li2p, Li3s, Li3d, Li4s, Li4p]
ListOfNames = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p']
for i in range(len(ListOfStates)):
    plt.plot(ListOfStates[i][0], ListOfStates[i][1], '-', label = ListOfNames[i])
plt.legend()
plt.title('Radial-wavefunctions for Li I, different states')
plt.xlabel('r [a.u]')
plt.ylabel(r'$P(r)$ $[a.u]$')
plt.xlim(0, 10)
plt.ylim(-0.5, 3)
plt.show()
print(Li1s[-1] * cm1converter)
print(Li2s[-1] * cm1converter)
print(Li2p[-1] * cm1converter)
print(Li3s[-1] * cm1converter)
print(Li3p[-1] * cm1converter)
print(Li3d[-1] * cm1converter)
print(Li4s[-1] * cm1converter)
print(Li4p[-1] * cm1converter)
"""
