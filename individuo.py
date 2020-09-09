from gene import Gene

class Individuo:
    F = 0.0
    CR = 0.0
    fitness = 0.0
    rmsd = 0.0
    genes = []
    pose = None
    index = 0
    distancias = []
    ss = []

    def __init__(self, F, CR):
        self.F = F
        self.CR = CR
        self.genes = []     

    def setF(self, f):
        self.F = f

    def getF(self):
        return self.F

    def setCR(self, cr):
        self.CR = cr

    def getCR(self):
        return self.CR

    def setFitness(self, fitness):
        self.fitness = fitness

    def getFitness(self):
        return self.fitness

    def setRmsd(self, rmsd):
        self.rmsd = rmsd
    
    def getRmsd(self):
        return self.rmsd

    def setPose(self, pose):
        self.pose = pose

    def getPose(self):
        return self.pose

    def setGenes(self, genes):
        self.genes = genes

    def getGenes(self):
        return self.genes

    def setIndex(self, index):
        self.index = index

    def getIndex(self):
        return self.index

    def setSS(self, ss):
        self.ss = ss

    def getSS(self):
        return self.ss