import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from numpy import polyfit
import math

global D1, D2, ub, limit, Convert

limit = 10 ** (-20)

Convert = 6.626 * 10 ** (-34) * 2.998 * 10 ** (8)/(1.602 * 10 ** (-19))

D2 = 358.05 #pixels

D1 = 221.073 #pixels

ub = 5.8 * 10 ** (-5) #eV/T

def Sigma(Da, Db):

    return 1/(2 * 3.085 * 10 ** (-3)) * ( Da ** 2 - Db ** 2) / (D2 ** 2 - D1 ** 2) * 2.998 * 10 ** 8 * 6.626 * 10 ** (-34)/(1.602*10**(-19))

def BField(Current):
    return Current * 0.1

CurrentList = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5] #Ampere
DaValues = [231.572, 244.285, 242.868, 253.953, 254.487, 265.238, 278.957, 280.347, 281.347, 283.003, 287.615, 314.590, 314.590, 314.590, 318.309, 322.040, 323.306, 331.110, 331.539, 342.539]
DbValues = [209.178, 193.743, 181.697, 178.473, 170.455, 161.709, 155.952, 143.231, 138.367, 147.544, 141.047, 127.422, 115.089, 114.300, 114.846, 133.873, 138.417, 138.922, 143.012, 152.539]
BfieldValues = np.array([BField(CurrentList[i]) for i in range(len(CurrentList))])
SigmaValues = np.array([Sigma(DaValues[i], DbValues[i]) for i in range(len(DaValues)) if len(DaValues) == len(DbValues)])

Model = linregress(BfieldValues, SigmaValues)
plt.plot(BfieldValues, Model[0] * BfieldValues + Model[1], label = 'Linear fit')
print(Model)

m_jg_j = Model[0] / ub
print(m_jg_j)

def PartialDa(Da, Db): #Da, Db are values
    partial = (Sigma(Da + limit, Db) - Sigma(Da, Db)) / limit
    return partial

def PartialDb(Da, Db): #Da, Db are values
    partial = (Sigma(Da, Db + limit) - Sigma(Da, Db))/ limit
    return partial

def PartialD1(Da, Db):
    partial = -D1/(3.085 * 10 ** (-3)) * (Da ** 2 - Db ** 2)/(D2 ** 2 - D1 ** 2)**2
    return partial

def PartialD2(Da, Db):
    partial = D2/(3.085 * 10 ** (-3)) * (Da ** 2 - Db ** 2)/(D2 ** 2 - D1 ** 2)**2
    return partial

ErrorValues = np.array([PartialDa(DaValues[i], DbValues[i])**2 * (DaValues[i] * 0.01)**2 + PartialDb(DaValues[i], DbValues[i])**2*(DbValues[i] * 0.01)**2 + PartialD1(DaValues[i], DbValues[i])**2 * (D1 * 0.01)**2
    + PartialD2(DaValues[i], DbValues[i])**2 * (D2 * 0.01)**2 for i in range(len(DaValues))])
RealErrorValues = (ErrorValues) ** (1/2) * Convert

plt.plot(BfieldValues, SigmaValues, '.', label = 'Data acquired')
plt.errorbar(BfieldValues, Model[0] * BfieldValues + Model[1], xerr = 0,  yerr = Model[4], linestyle = "None", label = "Error of Linear fit")
plt.errorbar(BfieldValues, SigmaValues, xerr = 0.1 * 0.005, yerr = RealErrorValues, linestyle = "None", label = "Error of measurement")
plt.legend()
plt.xlabel(r'Magnetic field $B$')
plt.ylabel(r'$\Delta \sigma$ [eV]')
plt.title(r"Transistion spectrum with Zeeman effect, with $\lambda = 470 nm$ filter")
plt.show()
