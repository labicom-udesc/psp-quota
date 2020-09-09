import random
import yaml
import bounds
from copy import copy
from log import Log
from individuo import Individuo
from gene import Gene
from jDE_services import jDEServices
from pathos.multiprocessing import ThreadingPool as Pool


class dsm_paralelo:
    
    def __init__(self, rosetta_conf):
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
                self.NP = config['NP']
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
        #self.NP = self.D*5
        self.C_min = (int)(self.NP ** 0.5)
        self.C_max = (int)(self.NP/10)
        self.target = None
        self.maxAval = self.NP * self.maxIteractions
        self.avaliacoes = 0
        self.bounds = bounds.Bounds()
        self.dssp = self.rosetta_conf.dssp
        self.caminho = 'testes_protocol/'+self.protocolo+'/'+self.proteinName+'/'
        self.nome = str(random.random())   # nome do arquivo para log


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
    def printPDB(self, pose):
        pose.dump_pdb(self.caminho+self.nome+".pdb")


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
        Log().best_score3(self.caminho+self.nome, self.pop[best_index].getFitness())
        Log().media_score3(self.caminho+self.nome, soma / self.NP)

    
    def fragment_insert(self, ind=None, n=100, temp=2.0, mode='3s'):

        pose = ind.getPose()
        #print('fit original %f \n'%ind.getFitness())
        
        mc = self.rosetta_conf.get_new_mc(pose, self.rosetta_conf.score3, temp)
        mover = self.rosetta_conf.get_mer(mode)
        trial = self.rosetta_conf.get_new_trial_mover(mover, mc)
        mover = trial        

        mc.set_temperature(temp)
        mover.apply(pose)

        pose.assign(mc.lowest_score_pose())
        #print('fit pose 1 %f \n'%mc.lowest_score())
        
        original = self.rosetta_conf.score3(pose)
        #print('fit original2 %f \n'%original)
        one_more = original is None
        evals = 0

        for _ in range(n):
            mover.apply(pose)
            #print('fit eval %i %f \n'%(evals,pose.energies().total_energy()))
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
        #print('ultimo fit %f \n'%ind.getFitness())
        ind = self.update_angle_from_pose(ind, pose)

        return ind

    
    def gera_trial(self, target, mut, f, cr):
        
        trial = Individuo(f, cr)
        pose = self.rosetta_conf.get_generalPose()
        genes = []
        ss_nova = ''
        ss = self.rosetta_conf.ssq3
        ss_target = target.getSS()
        ss_mut = mut.getSS()

        # loop para gerar novo indivíduo
        for j in range(0, self.D):
            gen = Gene()
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

            gen.setPhi(gen_doador.getPhi())
            gen.setPsi(gen_doador.getPsi())
            gen.setOmega(gen_doador.getOmega())
            gen.setSSType(gen_doador.getSSType())     
 
            # seta o objeto pose para calcular o fitness do individuo
            pose.set_omega((j+1), gen.getOmega())
            pose.set_phi((j+1), gen.getPhi())
            pose.set_psi((j+1), gen.getPsi())

            genes.append(gen)
            ss_nova += gen.getSSType()

        trial.setSS(ss_nova)            
        trial.setGenes(genes)
        trial.setPose(pose)
        trial.setFitness(self.rosetta_conf.score3(pose))
        self.avaliacoes += 1

        return trial


    def gera_trial2(self, seed, target):
        
        trial1 = copy(target)
        pose1 = target.getPose()
        trial2 = copy(seed)
        pose2 = seed.getPose()
        ss = self.rosetta_conf.ssq3        # SS prevista pelo preditor
        r = random.randint(0, self.D-9)

        # loop para gerar novo indivíduo
        for j in range(r, self.D):

            if ss[j] == ss[r]:

                # insere seed no target -> gera trial1
                trial1.genes[j].setPhi(seed.genes[j].getPhi())
                trial1.genes[j].setPsi(seed.genes[j].getPsi())
                trial1.genes[j].setOmega(seed.genes[j].getOmega())
                trial1.genes[j].setSSType(seed.genes[j].getSSType())
                pose1.set_phi((j+1), seed.genes[j].getPhi())
                pose1.set_psi((j+1), seed.genes[j].getPsi())
                pose1.set_omega((j+1), seed.genes[j].getOmega())

                # insere target no seed -> gera trial2
                trial2.genes[j].setPhi(target.genes[j].getPhi())
                trial2.genes[j].setPsi(target.genes[j].getPsi())
                trial2.genes[j].setOmega(target.genes[j].getOmega())
                trial2.genes[j].setSSType(target.genes[j].getSSType())
                pose2.set_phi((j+1), target.genes[j].getPhi())
                pose2.set_psi((j+1), target.genes[j].getPsi())
                pose2.set_omega((j+1), target.genes[j].getOmega())

            else:
                break
  
        trial1.setPose(pose1)
        trial2.setPose(pose2)
        trial1.setFitness(self.rosetta_conf.score3(pose1))
        trial2.setFitness(self.rosetta_conf.score3(pose2))
        self.avaliacoes += 2

        return trial1, trial2


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


    # Retorna 1 indice para um individou diferente de j e com especie tam i
    def get_ind_aleatorio(self, i, j):
        ind = random.randint(0, i-1)
        while ind == j or ind == 0:
            ind = random.randint(0, i-1)

        return ind 


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


    def geraMutacao(self, seed, seedA, ind, f, cr):
        Vmut = Individuo(f, cr)
        pose = self.rosetta_conf.get_generalPose()
        genes = []

        for j in range(0, self.D):
            gen = Gene()

            phi = seed.genes[j].getPhi() + (0.05 * (seedA.genes[j].getPhi() - ind.genes[j].getPhi()))
            psi = seed.genes[j].getPsi() + (0.05 * (seedA.genes[j].getPsi() - ind.genes[j].getPsi()))
            omega = seed.genes[j].getOmega() + (0.05 * (seedA.genes[j].getOmega() - ind.genes[j].getOmega()))
                                                   
            gen.setPhi(phi)
            gen.setPsi(psi)
            gen.setOmega(omega)

            pose.set_phi((j+1), phi)
            pose.set_psi((j+1), psi)
            pose.set_omega((j+1), omega)

            genes.append(gen)
                    
        Vmut.setGenes(genes)
                
        Vmut.setPose(pose)
        #Vmut.setFitness(self.rosetta_conf.score3(pose))
        self.dssp.apply(pose)    
        ss = pose.secstruct()
        Vmut.setSS(ss)
        self.avaliacoes += 1

        return Vmut


    # Retorna index de outra especie
    def get_outra_especie(self, index):
        r = random.randint(0, len(self.especies)-1)
        while self.especies[r][0].getIndex() == index:
            r = random.randint(0, len(self.especies)-1)
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

    
    def verifica_best_pop(self):
        best_index = 0
        bad_index = 0
        mediaFit = 0

        for i in range(0, self.NP):
            # Verifica o melhor indivíduo de toda a população 
            if self.pop[i].getFitness() <= self.pop[best_index].getFitness():
                best_index = i
            if self.pop[i].getFitness() >= self.pop[bad_index].getFitness():
                bad_index = i
            mediaFit = mediaFit + self.pop[i].getFitness()
            self.pop[i].setIndex(i)

        return best_index, bad_index, mediaFit
    

    def ciclo_evolutivo(self, especie):

            for j in range(0, len(especie)):  # colunas = individuos da especie
                seed = especie[0]                 # semente da espécie atual

                target = especie[j]      # individuo atual
                # Auto ajuste dos parâmetros jDE
                #F_new, CR_new = self.de_service.autoAjuste(self.target.getF(), self.target.getCR())
                F_new = self.de_service.ajuste_f()
                    
                # Retorna outra semente
                seedA = self.especies[self.get_outra_especie(seed.getIndex())][0]

                # Retorna outro individuo da mesma especie
                ind = especie[self.get_ind_aleatorio(len(especie), j)]

                ####### MUTAÇÃO 1 #######
                Vmut = self.geraMutacao(seed, seedA, ind, F_new, self.CR)

                ####### CROSSOVER #######
                trial = self.gera_trial(target, Vmut, F_new, self.CR)

                ####### MUTAÇÃO 2 #######
                if random.uniform(0,1) <= F_new:
                    trial = self.fragment_insert(ind=copy(trial), mode=self.get_mode())
                                        
                ####### SELEÇÃO #######
                if trial.getFitness() <= target.getFitness():
                    self.pop.append(copy(trial))
                else:
                    x1, x2 = self.conta_ss(trial, target)
                    if x1 > x2:
                        self.pop.append(copy(trial))
                    else:
                        self.pop.append(copy(target))


    def paraleliza(self, funcao: callable, i: int, pool_size: int):
        pool = Pool(4) # Pool de 4 threads

        pool.map(funcao, self.especies[i]) # divide population em 4 threads e executa o fitness para cada em paralelo
        #evaluated_population.sort(key=lambda x: x[1])

        #return np.array(evaluated_population, dtype=object)

        
    def otimiza(self):
        self.init_pop() # gera a população inicial 
                
        numGer = 1 # contador de gerações
        cont_best = 0 # contador para reinicialização das espécies (if reinit == 10)
        
        while numGer < self.maxIteractions and self.avaliacoes < self.maxAval:
            
            self.gera_vetor_dist()
            self.gera_especies()
            
            self.pop[:] = []  # limpa o vetor de população

            self.printEspecies()

            best_index = 0    # armazena o índice do melhor indivíduo da população geral
            bad_index = 0     # armazena o índice do pior indivíduo da população geral
            mediaFit = 0      # armazena o valor medio da população

            for i in range(0, len(self.especies)):         # linhas = especies
                self.paraleliza(self.ciclo_evolutivo, i, 4) # paraleliza os individuos dentro de cada especie
            
            best_index, bad_index, mediaFit = self.verifica_best_pop()

            Log().best_score3(self.caminho+self.nome, self.pop[best_index].getFitness())
            Log().media_score3(self.caminho+self.nome, mediaFit / self.NP)

            # conta ha quantas gerações o best não muda
            if self.pop[best_index].getFitness() >= self.min_fit:
                cont_best += 1
            else:
                cont_best = 0

            self.min_fit = self.pop[best_index].getFitness()      
            self.max_fit = self.pop[bad_index].getFitness()

            '''self.especies[:] = []
            self.especies = copy(especAux)
            especAux[:] = []'''

            #print('\nDIVERSIDADE %f \n' %self.avg_rmsd())
            numGer += 1

        print('---------------------------------SCORE3 %f' %self.pop[best_index].getFitness())
        self.printPDB(self.pop[best_index].getPose())

        return self.pop[best_index].getPose()