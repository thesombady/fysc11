def FineEnergy(n,j,z):
    Energy = -13.6 / n**2 * z**2 * (1 + (1/137*z)**2 / n ** 2 * (n /(j + 1/2) - 3/4))
    return Energy

EnergyDiff = FineEnergy(5, 1/2, 1) - FineEnergy(5, 3/2, 1)

print(EnergyDiff, 'eV')

def Convert(Value):
    Value = Value * 1.602 * 10**(-19)
    h = 6.626 * 10**(-34)
    return Value / h * 10**(-9)

print(Convert(EnergyDiff), 'GHz')
