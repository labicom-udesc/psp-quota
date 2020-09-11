import random
import yaml
import bounds
import time
from copy import copy
from log import Log
from individuo import Individuo
from gene import Gene
from jDE_services import jDEServices
from pathos.multiprocessing import ThreadingPool as Pool


class dsm2:
    
    def __init__(self, rosetta_conf, caminho, nome):
        with open("config.yaml", 'r') as stream:
            try:
                config = yaml.load(stream)
                #self.NP = config['NP']
                self.T1 = config['T1']
                self.T2 = config['T2']
                self.F_lower = config['F_lower']
                self.F_upper = config['F_upper']
                self.F = config['F']
                self.CR = config['CR']
                self.maxIteractions = config['maxIteractions']
                self.proteinName = config['proteinName']
                self.protocolo = config['protocolo']
                self.mutation = config['mutation']
                #self.NP = config['NP']
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
        self.target = None
        self.maxAval = self.NP * self.maxIteractions
        self.avaliacoes = 0
        self.bounds = bounds.Bounds()
        self.dssp = self.rosetta_conf.dssp
        #self.caminho = 'testes_protocol/'+self.protocolo+'/'+self.proteinName+'/'
        self.caminho = caminho
        #self.nome = str(random.random())   # nome do arquivo para log
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
        ind.setSS(ss)
        for k in range(0, self.D):
            gen = Gene()
            gen.setPhi(pose.phi(k + 1))
            gen.setPsi(pose.psi(k + 1))
            gen.setOmega(pose.omega(k + 1))
            gen.setSSType(ss[k])
            genes.append(gen)
        ind.setGenes(genes)

        return ind

    
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

        self.update_angle_from_pose(ind, pose)

        self.avaliacoes += 1
        ind.setPose(pose)
        return ind


    # Gera a população inicial
    def init_pop(self):
        best_index = 0
        bad_index = 0
        soma = 0

        for i in range(0, self.NP):
            ind = Individuo(self.F, self.CR)
            ind = self.new_indiv_frag(ind)     # gera indivíduos a partir dos fragmentos
            ind.setIndex(i)
            
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
        Log().best_score3(self.caminho+self.nome, self.avaliacoes, self.pop[best_index].getFitness())
        Log().media_score3(self.caminho+self.nome, self.avaliacoes, soma / self.NP)

    
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
        ind = self.update_angle_from_pose(ind, pose)

        return ind

    
    def gera_trial(self, target, mut, r):
        
        pose = target.getPose()
        ss = self.rosetta_conf.ssq3
        ss_target = target.getSS()
        ss_mut = mut.getSS()

        # loop para gerar novo indivíduo
        for j in range(r, self.D):

            if ss[j] == ss[r]:
                gen_doador = Gene()
                    
                if (ss_target[j] == ss[j] and ss_mut[j] == ss[j]) or (ss_target[j] != ss[j] and ss_mut[j] != ss[j]):
                    if random.uniform(0,1) <= 0.5:
                        gen_doador = target.genes[j]
                    else:
                        gen_doador = mut.genes[j]
                        gen_doador.setSSType(ss_mut[j])
                elif ss_target[j] == ss[j] and ss_mut[j] != ss[j]:
                    gen_doador = target.genes[j]
                else:
                    gen_doador = mut.genes[j]
                    gen_doador.setSSType(ss_mut[j])

                target.genes[j].setPhi(gen_doador.getPhi())
                target.genes[j].setPsi(gen_doador.getPsi())
                target.genes[j].setOmega(gen_doador.getOmega())

                pose.set_phi((j+1), gen_doador.getPhi())
                pose.set_psi((j+1), gen_doador.getPsi())
                pose.set_omega((j+1), gen_doador.getOmega())
            else:
                break                 

        target.setPose(pose)

        return target


    def gera_trial3(self, target, seed, seedA, ind):
        
        trial1 = copy(target)
        pose1 = target.getPose()
        ss = self.rosetta_conf.ssq3        # SS prevista pelo preditor
        r = random.randint(0, self.D-9)

        # loop para gerar novo indivíduo
        for j in range(r, self.D):

            if ss[j] == ss[r]:

                phi = seed.genes[j].getPhi() + (0.05 * (seedA.genes[j].getPhi() - ind.genes[j].getPhi()))
                psi = seed.genes[j].getPsi() + (0.05 * (seedA.genes[j].getPsi() - ind.genes[j].getPsi()))
                omega = seed.genes[j].getOmega() + (0.05 * (seedA.genes[j].getOmega() - ind.genes[j].getOmega()))
             
                # insere seed no target -> gera trial1
                trial1.genes[j].setPhi(phi)
                trial1.genes[j].setPsi(psi)
                trial1.genes[j].setOmega(omega)
                pose1.set_phi((j+1), phi)
                pose1.set_psi((j+1), psi)
                pose1.set_omega((j+1), omega)

            else:
                break
  
        trial1.setPose(pose1)
        trial1.setFitness(self.rosetta_conf.score3(pose1))
        self.dssp.apply(pose1)
        trial1.setSS(pose1.secstruct())
        self.avaliacoes += 1

        return trial1


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


    # Retorna 1 indivíduo aleatório da mesma especie, diferente de ind[i][j] e seed
    def get_ind_aleatorio(self, i, j):
        ind = random.randint(0, len(self.especies[i])-1)
        while ind == j or ind == 0:
            ind = random.randint(0, len(self.especies[i])-1)

        return self.especies[i][ind] 


    def forced_insertion_frag(self, ind):
        #score_before = ind.getFitness()
        #print('score antes %f \n' %score_before)
        pose = ind.getPose()

        self.rosetta_conf.get_mer(self.get_mode()).apply(pose)
        ind = self.update_angle_from_pose(ind, pose)
        ind.setFitness(self.rosetta_conf.score3(pose))

        self.avaliacoes += 1

        #score_after = ind.getFitness()
        #print('score depois %f \n' %score_after)

        return ind


    def geraMutacao(self, seed, seedA, ind, r):

        Vmut = copy(seed)
        pose = seed.getPose()
        ss = self.rosetta_conf.ssq3

        for j in range(r, self.D):

            if ss[j] == ss[r]:
                phi = seed.genes[j].getPhi() + (0.03 * (seedA.genes[j].getPhi() - ind.genes[j].getPhi()))
                psi = seed.genes[j].getPsi() + (0.03 * (seedA.genes[j].getPsi() - ind.genes[j].getPsi()))
                omega = seed.genes[j].getOmega() + (0.03 * (seedA.genes[j].getOmega() - ind.genes[j].getOmega()))
                                                   
                Vmut.genes[j].setPhi(phi)
                Vmut.genes[j].setPsi(psi)
                Vmut.genes[j].setOmega(omega)

                pose.set_phi((j+1), phi)
                pose.set_psi((j+1), psi)
                pose.set_omega((j+1), omega)
            else:
                break
                                 
        Vmut.setPose(pose)
        self.dssp.apply(pose)    
        ss = pose.secstruct()
        Vmut.setSS(ss)

        return Vmut


    # Retorna index de outra especie
    def get_outra_especie(self, tam, atual):
        r = random.randint(0, tam-1)
        while r == atual:
            r = random.randint(0, tam-1)
        return r


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

    def reinicia_populacao(self, ind):
        self.pop[:] = []
        for i in range(0, self.NP):
            novo = copy(ind)
            novo = self.forced_insertion_frag(novo)
            self.pop.append(novo)   

    def reinicia_especies(self):
        for i in range(0, len(self.especies)):
            tam = len(self.especies[i])
            corte = (int)(tam/3)
            while corte > 0:
                tam -= 1
                corte -=1
                ind = copy(self.especies[i][0])
                index = self.especies[i][tam].getIndex()
               # ind = Individuo(self.F, self.CR)
               # ind = self.new_indiv_frag(ind)
                ind = self.forced_insertion_frag(ind)
                ind.setIndex(index)
               # ind.setFitness(self.rosetta_conf.score3(ind.getPose()))
                self.especies[i][tam] = ind


    def conta_ss(self, trial, target):
        cont_trial = 0
        cont_target = 0
        ss = self.rosetta_conf.ssq3
        ss_trial = trial.getSS()
        ss_target = target.getSS()

        for j in range(0, self.D):
            if ss_trial[j] == ss[j]:
                cont_trial += 1
            if ss_target[j] == ss[j]:
                cont_target += 1

        return cont_trial, cont_target

        
    def otimiza(self):
        time_pop_inicial = time.time()

        self.init_pop() # gera a população inicial 
        self.gera_vetor_dist()
        self.gera_especies()
                
        numGer = 1 # contador de gerações

        time_ciclo_evolutivo_start = time.time()
        
        while numGer < self.maxIteractions and self.avaliacoes < self.maxAval:
                        
            self.pop[:] = []  # limpa o vetor de população
            #self.printEspecies()

            best_index = 0    # armazena o índice do melhor indivíduo da população geral
            bad_index = 0     # armazena o índice do pior indivíduo da população geral
            contInd = 0       # conta o numero de individuos na população
            mediaFit = 0      # armazena o valor medio da população

            for i in range(0, len(self.especies)):         # linhas = especies

                seed = self.especies[i][0]                 # semente da espécie atual
                
                for j in range(0, len(self.especies[i])):  # colunas = individuos da especie

                    self.target = self.especies[i][j]      # individuo atual

                    # Auto ajuste dos parâmetros jDE
                    F_new, CR_new = self.de_service.autoAjuste(self.target.getF(), self.target.getCR())
                                        
                    # Retorna outra semente
                    #seedA = self.especies[self.get_outra_especie(len(self.especies), i)][0]

                    # Retorna outro individuo da mesma especie
                    #ind = self.get_ind_aleatorio(i, j)
              
                    trial = copy(self.target)
                    trial.setF(F_new)
                    trial.setCR(CR_new)

                    ####### CROSSOVER #######
                    if random.uniform(0,1) <= CR_new:
                        r = random.randint(0, self.D-9)
                        #Vmut = self.geraMutacao(seed, seedA, ind, r)   ### MUTAÇÃO 1
                        trial = self.gera_trial(trial, seed, r)

                    ####### MUTAÇÃO 2 #######
                    trial = self.fragment_insert(ind=trial, mode=self.get_mode())
                                        
                    ####### SELEÇÃO #######
                    if trial.getFitness() < self.target.getFitness():
                        self.pop.append(copy(trial))
                    elif trial.getFitness() > self.target.getFitness():
                        self.pop.append(copy(self.target))
                    else:
                        x1, x2 = self.conta_ss(trial, self.target)
                        if x1 > x2:
                            self.pop.append(copy(trial))
                        else:
                            self.pop.append(copy(self.target))
                    self.pop[contInd].setIndex(contInd)     # seta o indice da pop geral

                    # Verifica o melhor indivíduo de toda a população 
                    if self.pop[contInd].getFitness() <= self.pop[best_index].getFitness():
                        best_index = contInd
                    if self.pop[contInd].getFitness() >= self.pop[bad_index].getFitness():
                        bad_index = contInd
                    mediaFit = mediaFit + self.pop[contInd].getFitness()

                    contInd += 1            
            
            Log().best_score3(self.caminho+self.nome, self.avaliacoes, self.pop[best_index].getFitness())
            Log().media_score3(self.caminho+self.nome, self.avaliacoes, mediaFit / self.NP)

            self.min_fit = self.pop[best_index].getFitness()      
            self.max_fit = self.pop[bad_index].getFitness()

            self.gera_vetor_dist()
            self.gera_especies()
            numGer += 1

        time_ciclo_evolutivo_end = time.time()

        Log().time(self.caminho+self.nome, 'Tempo ciclo evolutivo', time_ciclo_evolutivo_end - time_ciclo_evolutivo_start)
        Log().time(self.caminho+self.nome, 'Tempo pop inicial', time_ciclo_evolutivo_end - time_pop_inicial)

        rmsd_best = 0
        index_best = 0
        for s in range(0, len(self.especies)):
            nome = self.nome+'_'+str(s)
            self.printPDB(self.especies[s][0].getPose(), nome)
            rmsd = self.rosetta_conf.get_rmsd_from_native(self.especies[s][0].getPose())  
            Log().rmsd(self.caminho+self.nome, s, rmsd, self.especies[s][0].getFitness())
            #print('RMSD IND %i %f\n' %(s,rmsd))
            if rmsd_best == 0:
                rmsd_best = rmsd
                index_best = s
            elif rmsd < rmsd_best:
                rmsd_best = rmsd
                index_best = s  
        #print('melhor individuo %i %f\n' %(index_best, rmsd_best)) 

        return self.especies[index_best][0].getPose()