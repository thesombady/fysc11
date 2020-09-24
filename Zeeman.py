import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

global D1, D2, ub

D2 = 358.05 #um
D1 = 221.073 #um
ub = 5.8 * 10 ** (-5) #eV/T


def Sigma(Da, Db):
    return 1/(2 * 3.085) * ( Da ** 2 - Db ** 2) /(D2 ** 2 - D1 ** 2) * 0.00124

def BField(Current):
    return Current * 0.1

def LinRegress(Param):
    pass


CurrentList = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5]
DaValues = [231.572, 244.285, 242.868, 253.953, 254.487, 265.238, 278.957, 280.347, 281.347, 283.003, 287.615, 314.590, 314.590, 314.590, 318.309, 322.040, 323.306, 331.110, 331.539, 342.539]
DbValues = [209.178, 193.743, 181.697, 178.473, 170.455, 161.709, 155.952, 143.231, 138.367, 147.544, 141.047, 127.422, 115.089, 114.300, 114.846, 133.873, 138.417, 138.922, 143.012, 152.539]
BfieldValues = np.array([BField(CurrentList[i]) for i in range(len(CurrentList))])
SigmaValues = np.array([Sigma(DaValues[i], DbValues[i]) for i in range(len(DaValues)) if len(DaValues) == len(DbValues)])

XValues = BfieldValues.reshape((-1,1))

Model = LinRegress().fit(XValues, SigmaValues)
#r_sq = Model.score(x, y)
print(Model)
#print('coefficient of determination:', r_sq)


plt.plot(BfieldValues, SigmaValues, '.')
plt.xlabel(r'Magnetic field $B$')
plt.ylabel(r'$\Delta \sigma$ [eV]')
plt.title("Transistion spectrum with doppler effect")
#plt.show()
