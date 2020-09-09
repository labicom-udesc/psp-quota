import bounds
import random

class jDEServices:

    def __init__(self, T1, T2, F_lower, F_upper, rosetta_conf):
        self.T1 = T1
        self.T2 = T2
        self.F_lower = F_lower
        self.F_upper = F_upper
        self.rosetta_conf = rosetta_conf
        self.bounds = bounds.Bounds()
        self.q8_mufold = ['C', 'B', 'E', 'G', 'I', 'H', 'S', 'T']


    # Retorna a estrutura secundária para o aminoácido da posição k
    # de acordo com o preditor de estrutura secundária especificado
    def getSS_psipred(self, k):
        return self.rosetta_conf.psi_pred[k]


    def getSS(self, k):
        rand = random.uniform(0, 1)
        ss = self.rosetta_conf.ssq8[k]
        ss_prob = self.rosetta_conf.mufold_prob[k]
        x = 0.00
        for k in range(0, len(ss_prob)):
            x = x + float(ss_prob[k])
            if rand < x:
                return self.q8_mufold[k]

        return ss



    # j = dimensão/gene do indivíduo
    # angulo = phi=0, psi=1
    # limit = 'L', 'U'
    def getLimiteAngulo(self, j, angulo, limit):
        if self.getSS(j) == '':
            return 180
        elif limit == 'L':
            return self.bounds.getSecondaryLowerBound(self.getSS(j), angulo)
        elif limit == 'U':
            return self.bounds.getSecondaryUpperBound(self.getSS(j), angulo)


    def get_index(self, NP, i):
        r1 = random.randint(0, NP-1)
        while(r1 == i):
            r1 = random.randint(0, NP-1)

        r2 = random.randint(0, NP-1)
        while(r2 == i or r2 == r1):
            r2 = random.randint(0, NP-1)

        r3 = random.randint(0, NP-1)
        while(r3 == i or r3 == r2 or r3 == r1):
            r3 = random.randint(0, NP-1)

        return r1, r2, r3


    def autoAjuste(self, F, CR):
        if random.uniform(0, 1) < self.T1:
            F_new = self.F_lower + random.uniform(0,1)*self.F_upper
        else:
            F_new = F
        if random.uniform(0,1) < self.T2:
            CR_new = random.uniform(0,1)
        else:
            CR_new = CR
        return F_new, CR_new 


    def get_melhor_fitness(self, pop, r1, r2, r3):
        x1 = pop[r1].getFitness()
        x2 = pop[r2].getFitness()
        x3 = pop[r3].getFitness()
        if x1 < x2 and x1 < x3:
            return r1
        elif x2 < x1 and x2 < x3:
            return r2
        else:
            return r3


    def ajuste_cr(self, D, numGer, maxIteractions):
        return (40/D)**(-numGer/maxIteractions)

    def ajuste_f(self):
        return random.gauss(0.5, 0.3)



    


    