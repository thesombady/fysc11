
def FineEnergy(n,j,z):
    Energy = -13.6 / n**2 * z**2 * (1 + (1/137)**2 / n ** 2 * (n /(j + 1/2) - 3/4))
    return Energy

Energy = FineEnergy(5, 1, 48) - FineEnergy(5, 0, 48)

def WaveConvertion(Energy):
    WaveLenght = (Energy * 1.602*10**(-19) / (6.626 * 10 **(-34) * 2.99 * 10 **8)) ** (-1)
    return WaveLenght
print(WaveConvertion(Energy))
