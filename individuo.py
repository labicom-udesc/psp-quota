from gene import Gene
import pyrosetta
from copy import copy

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
    sum_ss = 0
    prob_ss = 0.0
    sum_cm = 0.0

    def __init__(self, F, CR):
        self.F = F
        self.CR = CR
        self.genes = []     
        self.pose = pyrosetta.Pose()

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
        self.pose.assign(pose)

    def getPose(self):
        return self.pose

    def setGenes(self, genes):
        self.genes = copy(genes)

    def getGenes(self):
        return self.genes

    def setIndex(self, index):
        self.index = index

    def getIndex(self):
        return self.index

    def setSS(self, ss):
        self.ss = copy(ss)

    def getSS(self):
        return self.ss

    def setSum_ss(self, sum_ss):
        self.sum_ss = sum_ss

    def getSum_ss(self):
        return self.sum_ss

    def setProb_ss(self, prob_ss):
        self.prob_ss = prob_ss

    def getProb_ss(self):
        return self.prob_ss

    def setSum_cm(self, sum_cm):
        self.sum_cm = sum_cm

    def getSum_cm(self):
        return self.sum_cm