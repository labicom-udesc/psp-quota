import random

class Gene:
    phi = 0.0
    psi = 0.0
    omega = 0.0
    r = []             # angulos da cadeia lateral
    name = ''          # nome do aminoacido = sigla de uma letra
    numAngles = 0      # numero de angulos da cadeira lateral
    ssType = ''        # letra que representa a estrutura secundária

    def __init__(self, x_min=-180.0, x_max=+180.0):
        self.r = []
        self.phi = x_min + random.uniform(0,1)*(x_max - x_min)
        self.psi = x_min + random.uniform(0,1)*(x_max - x_min)
        self.omega = 180.0

    def setPhi(self, phi):
        self.phi = phi

    def getPhi(self):
        return self.phi

    def setPsi(self, psi):
        self.psi = psi

    def getPsi(self):
        return self.psi

    def setOmega(self, omega):
        self.omega = omega

    def getOmega(self):
        return self.omega

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name

    def setSSType(self, ss):
        self.ssType = ss

    def getSSType(self):
        return self.ssType

    def setSideChain(self, r):
        for i in range(0, len(r)):
            self.r.append(r[i])

    def getSideChain(self):
        return self.r

    def setNumAngles(self, num):
        self.numAngles = num

    def getNumAngles(self):
        return self.numAngles