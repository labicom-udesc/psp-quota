import random
import yaml
import time
from copy import copy
from log import Log
from individuo import Individuo
from gene import Gene
from jDE_services import jDEServices
from repacker import Repacker
import pyrosetta
import statistics
import math


class dsm2:
    
    def __init__(self, rosetta_conf, caminho, nome):
        with open("config.yaml", 'r') as stream:
            try:
                config = yaml.load(stream)
                self.T1 = config['T1']
                self.T2 = config['T2']
                self.F_lower = config['F_lower']
                self.F_upper = config['F_upper']
                self.F = config['F']
                self.CR = config['CR']
                self.maxAval = config['maxAval']
                self.proteinName = config['proteinName']
                self.protocolo = config['protocolo']
            except yaml.YAMLError as exc:                
                print(exc)
        self.rosetta_conf = rosetta_conf
        self.pop = []
        self.especies = []
        self.min_fit = None
        self.max_fit = None
        self.frag_modes = ['3', '3s', '9', '9s']
        self.de_service = jDEServices(self.T1, self.T2, self.F_lower, self.F_upper, self.rosetta_conf)
        self.D = len(self.rosetta_conf.fasta) 
        self.NP = self.D*5
        self.C_min = (int)(self.NP ** 0.5)
        self.C_max = (int)(self.NP/10)
        self.c = 2
        self.m_nmdf = 0
        self.target = None
        self.avaliacoes = 0
        self.dssp = self.rosetta_conf.dssp
        self.caminho = caminho
        self.nome = nome


    # Gera os vetores de distância euclidiana
    # Cada individuo tem um vetor para armazenar as distancias
    ## Cada posição do vetor representa a posição do individuo dentro da população
    ### Logo, para o indivíduo i, a posição i do vetor terá valor 0
    def gera_vetor_dist(self):
        for i in range(0, len(self.pop)):
            distancias = []
            for j in range(0, len(self.pop)):
                if i==j:
                    distancias.append(0.0)
                elif j<i:
                    distancias.append(self.pop[j].distancias[i])
                else:
                    distancias.append(self.rosetta_conf.get_rmsd_from_pose(self.pop[i].getPose(), self.pop[j].getPose()))
            self.pop[i].distancias = distancias

    def avg_rmsd(self):
        s = 0
        c = 0

        for i in range(self.NP):
            for j in range(i + 1, self.NP):
                c += 1
                s += self.rosetta_conf.get_rmsd_from_pose(self.pop[i].getPose(), self.pop[j].getPose())

        return s / c

    
    def update_diversity(self):
        diversity = 0
        aux_1 = 0
        aux_2 = 0
        d = 0

        total_number_of_angles = 3 * self.D

        for i in range(0, self.NP):
            for j in range(i + 1, self.NP):
                aux_1 = 0

                ind_a = self.pop[i].genes
                ind_b = self.pop[j].genes

                for d in range(0, self.D):
                    aux_1 += (ind_a[d].phi - ind_b[d].phi) ** 2
                    aux_1 += (ind_a[d].psi - ind_b[d].psi) ** 2
                    aux_1 += (ind_a[d].omega - ind_b[d].omega) ** 2

                aux_1 = math.sqrt(aux_1) / total_number_of_angles

                if j == i + 1 or aux_2 > aux_1:
                    aux_2 = aux_1

            diversity += math.log(1.0 + aux_2)

        self.m_nmdf = max(self.m_nmdf, diversity)

        if self.m_nmdf > 0:
            diversity = diversity / self.m_nmdf
        else:
            diversity = 0.0

        return diversity


    # Retorna o individuo com o melhor fitness da população restante
    def get_ind_best_fit(self):
        best = 0
        for i in range(0, len(self.pop)):
            if self.pop[i].getFitness() <= self.pop[best].getFitness():
                best = i
        return self.pop[best]       
                

    # Imprime o arquivo pdb para um determinado Pose
    def printPDB(self, pose, nome):
        pose.dump_pdb(self.caminho+nome+".pdb")


    def update_angle_from_pose(self, ind, pose):
        genes = []
        self.dssp.apply(pose)    # populates the pose's Pose.secstruct
        ss = pose.secstruct()
        ss_pred = self.rosetta_conf.ss3   # SS predita pelo preditor
        ind.setSS(ss)
        sum_ss = 0
        prob_ss = 0.0
        d = 0
        for k in range(0, self.D):
            gen = Gene()
            gen.setPhi(pose.phi(k + 1))
            gen.setPsi(pose.psi(k + 1))
            gen.setOmega(pose.omega(k + 1))
            gen.setSSType(ss[k])
            genes.append(gen)

            if ss[k] == ss_pred[k]:
                sum_ss += 1

            if ss[k] == 'L':
                prob_ss = prob_ss + (float)(self.rosetta_conf.probs[d])
            elif ss[k] == 'H':
                prob_ss = prob_ss + (float)(self.rosetta_conf.probs[d+1])
            else:
                prob_ss = prob_ss + (float)(self.rosetta_conf.probs[d+2])
            d = d+3
        ind.setGenes(genes)
        ind.setSum_ss(sum_ss)
        ind.setProb_ss(prob_ss)

    
    def get_mode(self):
        # sorteia um modelo de fragmento a ser utilizado (3, 3s, 9, 9s)
        rand = random.randint(0, len(self.frag_modes)-1)
        mode = self.frag_modes[rand]
        return mode


    # Cria novo indivíduo com inserção de fragmentos tamanho '9' e score '0'
    def new_indiv_frag(self, ind, n=100, temp=2.0):
        pose = self.rosetta_conf.get_new_pose()        

        mc = self.rosetta_conf.get_new_mc(pose, self.rosetta_conf.score0, temp)

        seq = self.rosetta_conf.get_new_seq_mover()
        seq.add_mover(self.rosetta_conf.get_mer(self.frag_modes[2]))
        trial = self.rosetta_conf.get_new_trial_mover(seq, mc)
        folding = self.rosetta_conf.get_new_rep_mover(trial, n)

        folding.apply(pose)
        pose.assign(mc.lowest_score_pose())

        self.avaliacoes += 1
        ind.setPose(pose)


    # Gera a população inicial
    def init_pop(self):
        best_index = 0
        bad_index = 0
        soma = 0

        for i in range(0, self.NP):
            ind = Individuo(self.F, self.CR)
            ind.setIndex(i)
            self.new_indiv_frag(ind)     # gera indivíduos a partir dos fragmentos
            self.update_angle_from_pose(ind, ind.getPose())
                        
            ind.setFitness(self.rosetta_conf.score3(ind.getPose()))
            self.avaliacoes += 1
            self.pop.append(ind)

            # avalia indivíduos da população
            if self.pop[i].getFitness() <= self.pop[best_index].getFitness():
                best_index = i
            if self.pop[i].getFitness() >= self.pop[bad_index].getFitness():
                bad_index = i
            soma = soma + self.pop[i].getFitness()

        self.min_fit = self.pop[best_index].getFitness()      
        self.max_fit = self.pop[bad_index].getFitness()
        Log().best_score3(self.caminho+self.nome, self.pop[best_index].getFitness())

    
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

        self.avaliacoes += evals
        pose.assign(mc.lowest_score_pose())
        ind.setPose(pose)
        ind.setFitness(mc.lowest_score())
        self.update_angle_from_pose(ind, pose)

    
    def contact_map(self, target, trial):
        dist_trial = 0
        dist_target = 0
        d = 0
        poseTrial = trial.getPose()
        poseTarget = target.getPose()
        dTrial = 0.0
        dTarget = 0.0

        while d < len(self.rosetta_conf.map):
            res1 = (int)(self.rosetta_conf.map[d])
            res2 = (int)(self.rosetta_conf.map[d+1])
            prob = (float)(self.rosetta_conf.map[d+2])
            type_res1 = 'CB'
            type_res2 = 'CB'
            if self.rosetta_conf.fasta[res1-1] == 'G':
                type_res1 = 'CA'
            if self.rosetta_conf.fasta[res2-1] == 'G':
                type_res2 = 'CA'

            dTrial = poseTrial.residue(res1).xyz(type_res1).distance(poseTrial.residue(res2).xyz(type_res2))
            if dTrial <= 8.0:
                dist_trial = dist_trial + prob
            else:
                dist_trial = dist_trial + (prob/dTrial)

            dTarget = poseTarget.residue(res1).xyz(type_res1).distance(poseTarget.residue(res2).xyz(type_res2))
            if dTarget <= 8.0:
                dist_target = dist_target + prob
            else:
                dist_target = dist_target + (prob/dTarget)
            
            d = d+3
        target.setSum_cm(dist_target)
        trial.setSum_cm(dist_trial)


    def gera_trial2(self, trial, seed, r, CR_new):

        ss = self.rosetta_conf.ss3
        pose = trial.getPose()
        
        # loop para gerar novo indivíduo
        for j in range(r, self.D):
            if ss[j] == ss[r]:
                if random.uniform(0,1) <= CR_new:
                    trial.genes[j].setPhi(seed.genes[j].getPhi())
                    trial.genes[j].setPsi(seed.genes[j].getPsi())
                    trial.genes[j].setOmega(seed.genes[j].getOmega())
                    trial.genes[j].setSSType(seed.genes[j].getSSType())

                    pose.set_phi((j+1), seed.genes[j].getPhi())
                    pose.set_psi((j+1), seed.genes[j].getPsi())
                    pose.set_omega((j+1), seed.genes[j].getOmega())
            else:
                break

        trial.setPose(pose)


    def gera_trial(self, trial, seed, r, CR_new):

        ss = self.rosetta_conf.ss3
        pose = trial.getPose()
        
        # loop crossover uniforme
        for j in range(0, self.D):
            if random.uniform(0,1) <= CR_new:
                trial.genes[j].setPhi(seed.genes[j].getPhi())
                trial.genes[j].setPsi(seed.genes[j].getPsi())
                trial.genes[j].setOmega(seed.genes[j].getOmega())
                trial.genes[j].setSSType(seed.genes[j].getSSType())

                pose.set_phi((j+1), seed.genes[j].getPhi())
                pose.set_psi((j+1), seed.genes[j].getPsi())
                pose.set_omega((j+1), seed.genes[j].getOmega())

        trial.setPose(pose)


    def gera_trial3(self, trial, seed, r, CR_new):

        ss = self.rosetta_conf.ss3
        ss_trial = trial.getSS()
        ss_seed = seed.getSS()
        pose = trial.getPose()
        
        # loop crossover uniforme
        for j in range(0, self.D):
            if (ss_trial[j] == ss[j] and ss_seed[j] == ss[j]) or (ss_trial[j] != ss[j] and ss_seed[j] != ss[j]):
                if random.uniform(0,1) <= 0.5:
                    trial.genes[j].setPhi(seed.genes[j].getPhi())
                    trial.genes[j].setPsi(seed.genes[j].getPsi())
                    trial.genes[j].setOmega(seed.genes[j].getOmega())
                    trial.genes[j].setSSType(seed.genes[j].getSSType())

                    pose.set_phi((j+1), seed.genes[j].getPhi())
                    pose.set_psi((j+1), seed.genes[j].getPsi())
                    pose.set_omega((j+1), seed.genes[j].getOmega())
            
            else:
                if ss_seed[j] == ss[j]:
                    trial.genes[j].setPhi(seed.genes[j].getPhi())
                    trial.genes[j].setPsi(seed.genes[j].getPsi())
                    trial.genes[j].setOmega(seed.genes[j].getOmega())
                    trial.genes[j].setSSType(seed.genes[j].getSSType())

                    pose.set_phi((j+1), seed.genes[j].getPhi())
                    pose.set_psi((j+1), seed.genes[j].getPsi())
                    pose.set_omega((j+1), seed.genes[j].getOmega())

        trial.setPose(pose)

    # Procura na população restante o indivíduo com menor distância da semente atual
    def get_ind_menor_dist(self, seed):
        menor_index = self.pop[0].getIndex()
        ind = self.pop[0]
        for i in range(1, len(self.pop)):
            if seed.distancias[self.pop[i].getIndex()] < seed.distancias[menor_index]:
                menor_index = self.pop[i].getIndex()
                ind = self.pop[i]
        return ind


    def get_tam_especie(self, seed):
        if (int)(self.min_fit) == (int)(self.max_fit):
            s = (int)(self.NP/10)
        else:
            s = (int)((((seed.getFitness()-self.max_fit)/(self.min_fit - self.max_fit))*(self.C_max - self.C_min)) + self.C_min)
        return s


    # Divide a população em espécies e armazena numa matriz, onde cada linha eh uma espécie
    # Ao final desse processo, não haverá mais indivíduos na matriz população
    def gera_especies(self):
        self.especies[:] = []
        cont = self.NP
        while cont > 0:
            especie = []
            seed = self.get_ind_best_fit()    # semente
            especie.append(seed)              # adiciona a semente na especie
            self.pop.remove(seed)             # remove a semente da população
            cont -= 1
            s = self.get_tam_especie(seed)
            if cont-s < self.C_min:
                s = cont + 1
            for i in range(0, s-1):
                ind = self.get_ind_menor_dist(seed)
                especie.append(ind)
                self.pop.remove(ind)   
                cont -= 1         
            self.especies.append(especie)


    def printPop(self):
        print('POPULACAO \n')
        for i in range(0,self.NP):
            print('ind %i fit %f\n' %(self.pop[i].getIndex(), self.pop[i].getFitness()))


    def printDistancias(self):
        print('DISTANCIAS\n')
        for i in range(0,self.NP):
            print(self.pop[i].distancias)
            print('\n')


    def printEspecies(self):
        print('ESPECIES\n')
        for i in range(0, len(self.especies)):
            print('especie %i \n' %i)
            for j in range(0, len(self.especies[i])):
                if j==0:
                    print(self.especies[i][j].getSS())
                print('%i %f' %(self.especies[i][j].getIndex(), self.especies[i][j].getFitness()))
            break


    def calcula_psc(self):
        result = 0
        Ess_pop = []
        Ess_mean = 0.0
        Ess_var = 0.0

        for i in range(0, self.NP):
            Ess_pop.append(self.pop[i].getSum_ss())

        Ess_mean = statistics.mean(Ess_pop)
        Ess_var = statistics.pvariance(Ess_pop, Ess_mean)
        val = (-self.c)*Ess_var*(math.pow(((Ess_mean/self.D)-1),2))
        result = math.exp(val)

        #print('media ESS %f %f %f\n' %(Ess_mean, Ess_var, result))

        return result

        
    def otimiza(self):

        time_start = time.time()

        self.init_pop() # gera a população inicial 
                
        numGer = 1 # contador de gerações
        div_usado = 0
        
        while self.avaliacoes < self.maxAval:
                        
            psc = self.calcula_psc()
            self.gera_vetor_dist()
            self.gera_especies()

            self.pop[:] = []  # limpa o vetor de população

            best_index = 0    # armazena o índice do melhor indivíduo da população geral
            bad_index = 0     # armazena o índice do pior indivíduo da população geral
            contInd = 0       # conta o numero de individuos na população

            for i in range(0, len(self.especies)):         # linhas = especies

                seed = self.especies[i][0]                 # semente da espécie atual
                
                for j in range(0, len(self.especies[i])):  # colunas = individuos da especie

                    self.target = self.especies[i][j]      # individuo atual

                    # Auto ajuste dos parâmetros jDE
                    F_new, CR_new = self.de_service.autoAjuste(self.target.getF(), self.target.getCR())
                                                      
                    trial = Individuo(F_new, CR_new)
                    trial.setGenes(self.target.getGenes())
                    trial.setPose(self.target.getPose())
                    trial.setSS(self.target.getSS())

                    ####### CROSSOVER #######
                    r = random.randint(0, self.D-3)
                    self.gera_trial2(trial, seed, r, CR_new)

                    ####### MUTAÇÃO #######
                    self.fragment_insert(ind=trial, mode=self.get_mode())
                    
                    ####### SELEÇÃO #######

                    if random.uniform(0,1) <= psc:
                        # seleção baseada em contato
                        self.contact_map(self.target, trial)
                        if trial.getSum_cm() > self.target.getSum_cm():
                            self.pop.append(trial)
                        else:
                            self.pop.append(self.target)
                    else:
                        # seleção baseada em estrutura
                        if trial.getProb_ss() > self.target.getProb_ss():
                            self.pop.append(trial)
                        else:
                            self.pop.append(self.target)

                    self.pop[contInd].setIndex(contInd)     # seta o indice da pop geral

                    # Verifica o melhor indivíduo de toda a população 
                    if self.pop[contInd].getFitness() <= self.pop[best_index].getFitness():
                        best_index = contInd
                    if self.pop[contInd].getFitness() >= self.pop[bad_index].getFitness():
                        bad_index = contInd

                    contInd += 1            
            
            # Log: registra energia a cada 10mil avaliações
            div = (int)(self.avaliacoes/10000)
            if div > div_usado:
                #media_rmsd = self.avg_rmsd()
                #diversidade = self.update_diversity()
                #Log().rmsd_diversity(self.caminho+self.nome, media_rmsd, diversidade, psc)
                Log().best_score3(self.caminho+self.nome, self.pop[best_index].getFitness())
                div_usado = div

            self.min_fit = self.pop[best_index].getFitness()      
            self.max_fit = self.pop[bad_index].getFitness()

            numGer += 1
        
        rmsd_best = 0
        index_best = 0

        repack = Repacker(self.rosetta_conf)

        for s in range(0, len(self.pop)):
            nome = self.nome+'_'+str(s)
            rmsd = self.rosetta_conf.get_rmsd_from_native(self.pop[s].getPose())
            if rmsd_best == 0:
                rmsd_best = rmsd
                index_best = s
            elif rmsd < rmsd_best:
                rmsd_best = rmsd
                index_best = s

        
        pose = pyrosetta.Pose()
        
        pose.assign(self.pop[index_best].getPose())
        pose = repack.repack(pose)
    
        score  = self.rosetta_conf.scorefxn(pose)
        rmsd = self.rosetta_conf.get_rmsd_from_native(pose) 
        
        gdt = self.rosetta_conf.get_gdt(pose)

        time_end = time.time()

        Log().repack(self.caminho+self.proteinName, self.nome, rmsd, gdt, score, time_end - time_start)
        
        print('gdt %f\n' %gdt)
        pose.dump_pdb(self.caminho+self.nome+"_repacked.pdb")