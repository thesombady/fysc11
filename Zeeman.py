import matplotlib.pyplot as plt
import numpy as np


"""
def FineEnergy(n,j,z):
    Energy = -13.6 / n**2 * z**2 * (1 + (1/137)**2 / n ** 2 * (n /(j + 1/2) - 3/4))
    return Energy

def WaveConvertion(Energy):
    WaveLenght = (Energy * 1.602*10**(-19) / (6.626 * 10 **(-34) * 2.99 * 10 **8)) ** (-1)
    return WaveLenght

Energy = FineEnergy(5, 2, 2) - FineEnergy(6, 1, 2)
print(WaveConvertion(Energy))
print(Energy)
#Energy = FineEnergy(5, 1/2, 2) - FineEnergy(6, 1/2, 2)
#print(WaveConvertion(Energy))
#print(Energy)
"""
global D1, D2, ub

D1 = 0.5309 #um
D2 = 0.8827 #um
ub = 5.8 * 10 ** (-5) #eV/T
#D1 = 215 #pixel
#D2 = 352 #pixel


def Sigma(Da, Db):
    return 1/(2 * 3.085) * ( Da ** 2 - Db ** 2) /(D1 ** 2 - D2 ** 2)

def BField(Current):
    return Current * 0.1

CurrentList = [0.25, 0.5, 0.75, 1, 1.25 , 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75]#, 4, 4.25, 4.5, 4,75, 5]
#           0.25    0.5     0.75    1       1.25    1.5     1.75    2       2.25    2.5     2.75    3       3.25    3.5     3.75
DaValues = [0.9132, 0.5962, 0.6020, 0.5881, 0.6403, 0.6586, 0.6689, 0.6696, 0.6912, 0.8058, 0.7761, 0.7557, 0.7466, 0.7601, 0.7501]
DbValues = [0.5554, 0.5020, 0.4966, 0.4629, 0.4642, 0.4301, 0.4149, 0.3911, 0.3684, 0.7303, 0.7136, 0.5738, 0.5572, 0.5831, 0.5798]

BfieldValues = [BField(CurrentList[i]) for i in range(len(CurrentList))]
SigmaValues = [Sigma(DbValues[i], DaValues[i]) for i in range(len(DaValues)) if len(DaValues) == len(DbValues)]

plt.plot(BfieldValues, SigmaValues, '.')
plt.show()
