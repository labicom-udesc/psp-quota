import random
import yaml
import bounds
from copy import copy
from log import Log
from individuo import Individuo
from gene import Gene
from jDE_services import jDEServices


class jDE:
    
    def __init__(self, rosetta_conf):
        with open("config.yaml", 'r') as stream:
            try:
                config = yaml.load(stream)
                self.NP = config['NP']
                self.T1 = config['T1']
                self.T2 = config['T2']
                self.F_lower = config['F_lower']
                self.F_upper = config['F_upper']
                self.F = config['F']
                self.CR = config['CR']
                self.maxIteractions = config['maxIteractions']
                self.proteinName = config['proteinName']
                self.protocolo = config['protocolo']
            except yaml.YAMLError as exc:                
                print(exc)
        self.rosetta_conf = rosetta_conf
        self.pop = []
        self.target = None
        self.mutation = None
        self.trial = None
        self.avaliacoes = self.NP * self.maxIteractions
        self.frag_modes = ['3', '3s', '9', '9s']
        self.de_service = jDEServices(self.T1, self.T2, self.F_lower, self.F_upper, self.rosetta_conf)
        self.D = len(self.rosetta_conf.fasta) 
        self.bounds = bounds.Bounds()
        self.caminho = 'testes_protocol/'+self.protocolo+'/'+self.proteinName+'/'
        self.nome = str(random.random())   # nome do arquivo para log


    # Imprime o arquivo pdb para um determinado Pose
    def printPDB(self, pose):
        pose.dump_pdb(self.caminho+self.nome+".pdb")


    def new_indiv_SS(self, ind): 
        cont = 0       
        pose = self.rosetta_conf.get_new_pose()
        while cont < 100:
            score = 0            
            for k in range(0, self.D):
                gen = Gene()
                ss = self.de_service.getSS(k)
                gen.setSSType(ss)
                gen.setName(self.rosetta_conf.fasta[k])
                if ss != '':
                    phi, psi, omega = self.bounds.generateRandomAngles(ss)
                    gen.setPhi(phi)
                    gen.setPsi(psi)
                    gen.setOmega(omega)

                pose.set_phi(k + 1, gen.getPhi())
                pose.set_psi(k + 1, gen.getPsi())
                pose.set_omega(k + 1, gen.getOmega())

                ind.genes.append(gen)

            score = self.rosetta_conf.score0(pose)
            cont = cont + 1
            if int(score) == 0:
                break
        
        ind.setPose(pose)
        return ind   


    def update_angle_from_pose(self, ind, pose):
        genes = []
        for k in range(0, self.D):
            gen = Gene()
            gen.setPhi(pose.phi(k + 1))
            gen.setPsi(pose.psi(k + 1))
            gen.setOmega(pose.omega(k + 1))
            genes.append(gen)
        ind.setGenes(genes)

        return ind

    
    def get_mode(self):
        # sorteia um modelo de fragmento a ser utilizado (3, 9, 3s, 9s)
        rand = random.randint(0, len(self.frag_modes)-1)
        mode = self.frag_modes[rand]
        return mode


    # Cria novo indivíduo com inserção de fragmentos
    def new_indiv_frag(self, ind, n=100, temp=2.0):
        pose = self.rosetta_conf.get_new_pose()        

        mc = self.rosetta_conf.get_new_mc(pose, self.rosetta_conf.score0, temp)

        seq = self.rosetta_conf.get_new_seq_mover()
        seq.add_mover(self.rosetta_conf.get_mer(self.get_mode()))
        trial = self.rosetta_conf.get_new_trial_mover(seq, mc)
        folding = self.rosetta_conf.get_new_rep_mover(trial, n)

        folding.apply(pose)
        pose.assign(mc.lowest_score_pose())

        self.update_angle_from_pose(ind, pose)

        self.avaliacoes -= 1
        ind.setPose(pose)
        return ind


    # Gera a população inicial
    def init_pop(self):
        best_index = 0
        soma = 0

        for i in range(0, self.NP):
            ind = Individuo(self.F, self.CR)
            #ind = self.new_indiv_frag(ind)
            ind = self.new_indiv_SS(ind)
            
            ind.setFitness(self.rosetta_conf.score3(ind.getPose()))
            self.pop.append(ind)

            # avalia indivíduos da população
            if self.pop[i].getFitness() <= self.pop[best_index].getFitness():
                best_index = i
            soma = soma + self.pop[i].getFitness()

        Log().best_score3(self.caminho+self.nome, self.pop[best_index].getFitness())
        Log().media_score3(self.caminho+self.nome, soma / self.NP)


    def geraMutacao(self, r1, r2, r3, j):
        # pega individuos aleatorios para fazer a mutação

        phi = self.pop[r1].genes[j].getPhi() + (self.trial.getF() * (self.pop[r2].genes[j].getPhi() - self.pop[r3].genes[j].getPhi()))
        psi = self.pop[r1].genes[j].getPsi() + (self.trial.getF() * (self.pop[r2].genes[j].getPsi() - self.pop[r3].genes[j].getPsi()))
        omega = self.pop[r1].genes[j].getOmega() + (self.trial.getF() * (self.pop[r2].genes[j].getOmega() - self.pop[r3].genes[j].getOmega()))
                                       
        # check bounds
        if phi < self.de_service.getLimiteAngulo(j, 0, 'L'):
            phi = (self.de_service.getLimiteAngulo(j, 0, 'L') + phi) / 2.0
        elif phi > self.de_service.getLimiteAngulo(j, 0, 'U'):
            phi = (self.de_service.getLimiteAngulo(j, 0, 'U') + phi) / 2.0
                    
        if psi < self.de_service.getLimiteAngulo(j, 1, 'L'):
            psi = (self.de_service.getLimiteAngulo(j, 1, 'L') + psi) / 2.0
        elif psi > self.de_service.getLimiteAngulo(j, 1, 'U'):
            psi = (self.de_service.getLimiteAngulo(j, 1, 'U') + psi) / 2.0
                    
        if omega > 180:
            omega -= 360

        return phi, psi, omega

    
    def fragment_insert(self, ind=None, n=100, temp=2.0, mode='3s'):

        pose = ind.getPose()
        
        mc = self.rosetta_conf.get_new_mc(pose, self.rosetta_conf.score3, temp)
        mover = self.rosetta_conf.get_mer(mode)
        trial = self.rosetta_conf.get_new_trial_mover(mover, mc)
        mover = trial        

        mc.set_temperature(temp)
        mover.apply(pose)

        pose.assign(mc.lowest_score_pose())
        
        original = self.rosetta_conf.score3(pose)
        one_more = original is None
        evals = 0

        for _ in range(n):
            mover.apply(pose)
            evals += 1
            if pose.energies().total_energy() < original:
                if not one_more:
                    break
                else:
                    one_more = False

        self.avaliacoes -= evals
        pose.assign(mc.lowest_score_pose())
        ind.setPose(pose)
        ind.setFitness(mc.lowest_score())
        ind = self.update_angle_from_pose(ind, pose)

        return ind


    def gera_trial(self, i):
        
        pose = self.rosetta_conf.get_generalPose()

        #r1, r2, r3 = self.de_service.get_index(self.NP, i)
        genes = []

        # loop para gerar novo indivíduo
        for j in range(0, self.D):
            r = random.uniform(0,1)
            gen = Gene()
                    
            if r <= self.trial.getCR() or j == random.randint(0, self.D):
                phi, psi, omega = self.mutation.genes[j].getPhi(), self.mutation.genes[j].getPsi(), self.mutation.genes[j].getOmega()
                gen.setPhi(phi)
                gen.setPsi(psi)
                gen.setOmega(omega)
            else:
                gen.setPhi(self.target.genes[j].getPhi())
                gen.setPsi(self.target.genes[j].getPsi())
                gen.setOmega(self.target.genes[j].getOmega())
 
            # seta o objeto pose para calcular o fitness do individuo
            pose.set_omega((j+1), gen.getOmega())
            pose.set_phi((j+1), gen.getPhi())
            pose.set_psi((j+1), gen.getPsi())

            genes.append(gen)
                    
        self.trial.setGenes(genes)
                
        self.trial.setPose(pose)
        self.trial.setFitness(self.rosetta_conf.score3(pose))
        self.avaliacoes -= 1

    
    def otimiza(self):
        popAux = []

        self.init_pop() # gera a população inicial

        m = 0
                    
        while m < self.maxIteractions and self.avaliacoes > 0:

            best_index = 0
            soma = 0
            
            for i in range(0, self.NP):
                self.target = self.pop[i]   # pega o indivíduo atual

                F_new, CR_new = self.de_service.autoAjuste(self.target.getF(), self.target.getCR())

                # cria um vetor de mutação com base no individuo target
                self.mutation = self.fragment_insert(ind=copy(self.target), mode=self.get_mode())

                self.trial = Individuo(self.F, self.CR)                
                self.trial.setF(F_new)
                self.trial.setCR(CR_new)

                self.gera_trial(i)

                # verifica se o individuo trial é melhor que target e add na nova população                
                if self.trial.getFitness() < self.target.getFitness():
                    popAux.append(self.trial)
                else:
                    popAux.append(self.target)
                
                # avalia os indivíduos da população 
                if popAux[i].getFitness() <= popAux[best_index].getFitness():
                    best_index = i
                soma = soma + popAux[i].getFitness()
            
            # salva log dos indivíduos a cada 10 gerações
            if m % 10 == 0:
                Log().best_score3(self.caminho+self.nome, popAux[best_index].getFitness())
                Log().media_score3(self.caminho+self.nome, soma / self.NP)

            self.pop[:] = []
            self.pop = copy(popAux)
            popAux[:] = []
            m += 1

        print('---------------------------------SCORE3 %f' %self.pop[best_index].getFitness())
        self.printPDB(self.pop[best_index].getPose())

        return self.pop[best_index].getPose()
