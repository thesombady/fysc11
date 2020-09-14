import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

def BindingEnergy(A, Z, Pair = 1):
    c1 = 15.9
    c2 = 18.4
    c3 = 0.71
    c4 = 23.2
    c5 = 11.5
    return c1 * A - c2 * A **(2/3) - c3 * Z**2 / (A**(1/3)) - c4 * (A - 2*Z)**2 / A + c5 * Pair/ (A)**(1/2)
for i in range(0, 44, 2): #Calcium 62 does not exits because it's below neutron dripline.
    print(BindingEnergy(i + 20, 20) - BindingEnergy(i + 19, 20), i + 20)
print("---------------------------------------------")
for i in range(0, 44, 2):#36
    print(BindingEnergy(i + 20, 20) - BindingEnergy(i + 19, 20 - 1), i + 20)

def ColumbEnergy(Z):
    K = 1
    return 3/5 * K * ((Z**2)-(Z-1)**2)**(2/3)
print("---------------------------------------------")
print(BindingEnergy(11, 5, 0) - BindingEnergy(11, 6, 0))

EnergyValue = np.array([0.97, 1.2, 1.73, 1.74, 3.03, 3.25, 3.92, 4.38, 4.94, 5.51, 5.47])

Nucleinumber = [(6, 11), (7, 13), (8, 15), (9, 17), (12, 23), (13, 25), (15,29), (16, 31), (18,35), (20,39), (21, 41)]
"""
XAxis = []
for i in range(len(Nucleinumber)):
    Value = Nucleinumber[i][1]
    Energy = ColumbEnergy(Value)
    XAxis.append(Energy)
"""
XAxis = np.array([ColumbEnergy(Nucleinumber[i][1]) for i in range(len(Nucleinumber))])



Regression = stats.linregress(XAxis, EnergyValue)
print(1/Regression[0])
print(Regression)
plt.title('Beta - decay')
plt.plot(XAxis, EnergyValue, marker = '.', linestyle = '', label = 'Orignal Data')
plt.plot(XAxis, Regression[1] + Regression[0] * XAxis, label = 'Fitted')
plt.text(6, 1, f'Error of fitted line is {Regression[4]}')
plt.text(7, 2, f'y = {1/Regression[0]}x + {Regression[1]}')
plt.legend(loc = 'best')
plt.xlabel('Nuclei')
plt.ylabel('Q(EC)')
plt.show()
