import random
import yaml
import bounds
from copy import copy
from log import Log
from individuo import Individuo
from gene import Gene
from jDE_services import jDEServices


class dsmDE:
    
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
            except yaml.YAMLError as exc:                
                print(exc)
        self.rosetta_conf = rosetta_conf
        self.pop = []
        self.especies = []
        self.min_fit = None
        self.max_fit = None        
        self.n_aval = 0
        self.frag_modes = ['3', '3s', '9', '9s']
        self.de_service = jDEServices(self.T1, self.T2, self.F_lower, self.F_upper, self.rosetta_conf)
        self.D = len(self.rosetta_conf.fasta) 
        self.NP = self.D*5
        self.C_min = (int)(self.NP ** 0.5)
        self.C_max = (int)(self.NP/10)
        self.target = None
        self.avaliacoes = self.NP * self.maxIteractions
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
                    #distancias.append(self.dist_euclidiana(self.pop[i], self.pop[j]))
                    #print('metodo %f \n'%self.dist_euclidiana(self.pop[i], self.pop[j]))
                    distancias.append(self.rosetta_conf.get_rmsd_from_pose(self.pop[i].getPose(), self.pop[j].getPose()))
                    #print('rmsd %f \n'%self.rosetta_conf.get_rmsd_from_pose(self.pop[i].getPose(), self.pop[j].getPose()))
            self.pop[i].distancias = distancias


    # Calcula a distancia euclidiana entre 2 indivíduos
    def dist_euclidiana(self, ind1, ind2):
        soma = 0.00
        for d in range(0, self.D):
            phi = (ind1.genes[d].getPhi() - ind2.genes[d].getPhi()) ** 2
            psi = (ind1.genes[d].getPsi() - ind2.genes[d].getPsi()) ** 2
            omega = (ind1.genes[d].getOmega() - ind2.genes[d].getOmega()) ** 2
            soma = soma + phi + psi + omega
        soma = soma ** 0.5
        return soma


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
        self.dssp.apply(pose)    # populates the pose's Pose.secstruct
        ss = pose.secstruct()
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
        # sorteia um modelo de fragmento a ser utilizado (3, 9, 3s, 9s)
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

        self.avaliacoes -= 1
        self.n_aval += 1
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
            #ind = self.new_indiv_SS(ind)
            ind.setIndex(i)
            
            ind.setFitness(self.rosetta_conf.score3(ind.getPose()))
            self.avaliacoes -= 1
            self.n_aval += 1
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


    def check_bounds(self, phi, psi, omega, j):
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


    # Gera mutação apenas na dimensão j, segundo o parâmetro de cruzamento CR
    def gera_mutacao(self, seed, F_new, seedA, indA, j):

        phi = seed.genes[j].getPhi() + (F_new * (seedA.genes[j].getPhi() - indA.genes[j].getPhi()))
        psi = seed.genes[j].getPsi() + (F_new * (seedA.genes[j].getPsi() - indA.genes[j].getPsi()))
        omega = seed.genes[j].getOmega() + (F_new * (seedA.genes[j].getOmega() - indA.genes[j].getOmega()))
                                       
        # check bounds
        phi, psi, omega = self.check_bounds(phi, psi, omega, j)  

        return phi, psi, omega

    
    def fragment_insert(self, ind=None, n=25, temp=2.0, mode='3s'):

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

        self.avaliacoes -= evals
        self.n_aval += evals
        pose.assign(mc.lowest_score_pose())
        ind.setPose(pose)
        ind.setFitness(mc.lowest_score())
        #print('ultimo fit %f \n'%ind.getFitness())
        ind = self.update_angle_from_pose(ind, pose)

        return ind

    def gera_trial_2(self, seed, ind1, ind2, F_new, CR_new):
        trial = Individuo(F_new, CR_new)
        pose = self.target.getPose()
        genes = self.target.getGenes()
        troca = []

        # Seleciona os indices para troca
        troca.append(random.randint(0, self.D-3))
        troca.append(random.randint(0, self.D-3))
        while troca[1] == troca[0]:
            troca[1] = random.randint(0, self.D-3)
        troca.append(random.randint(0, self.D-3))
        while troca[2] == troca[0] or troca[2] == troca[0]:
            troca[2] = random.randint(0, self.D-3)

        # Mutação seed
        for j in range(0,2):    
            gen = Gene()
            gen.setPhi(seed.genes[troca[0]+j].getPhi())
            gen.setPsi(seed.genes[troca[0]+j].getPsi())
            gen.setOmega(seed.genes[troca[0]+j].getOmega())  
            pose.set_omega((troca[0]+j+1), gen.getOmega())
            pose.set_phi((troca[0]+j+1), gen.getPhi())
            pose.set_psi((troca[0]+j+1), gen.getPsi())
            genes[troca[0]+j] = gen

        # Mutação ind1
        for j in range(0,2):    
            gen = Gene()
            gen.setPhi(ind1.genes[troca[1]+j].getPhi())
            gen.setPsi(ind1.genes[troca[1]+j].getPsi())
            gen.setOmega(ind1.genes[troca[1]+j].getOmega())  
            pose.set_omega((troca[1]+j+1), gen.getOmega())
            pose.set_phi((troca[1]+j+1), gen.getPhi())
            pose.set_psi((troca[1]+j+1), gen.getPsi())
            genes[troca[1]+j] = gen

        # Mutação ind2
        for j in range(0,2):    
            gen = Gene()
            gen.setPhi(ind2.genes[troca[2]+j].getPhi())
            gen.setPsi(ind2.genes[troca[2]+j].getPsi())
            gen.setOmega(ind2.genes[troca[2]+j].getOmega())  
            pose.set_omega((troca[2]+j+1), gen.getOmega())
            pose.set_phi((troca[2]+j+1), gen.getPhi())
            pose.set_psi((troca[2]+j+1), gen.getPsi())
            genes[troca[2]+j] = gen

        trial.setGenes(genes)         
        trial.setPose(pose)
        trial.setFitness(self.rosetta_conf.score3(pose))
        self.avaliacoes -= 1
        self.n_aval += 1
        return trial


    def gera_trial_3(self, seed, F_new, CR_new):
        trial = Individuo(F_new, CR_new)
        pose = seed.getPose()
        genes = seed.getGenes()

        # Seleciona os indices para troca
        r = random.randint(0, self.D-3)

        # Mutação seed
        for j in range(0,2):    
            gen = Gene()
            gen.setPhi(self.target.genes[r+j].getPhi())
            gen.setPsi(self.target.genes[r+j].getPsi())
            gen.setOmega(self.target.genes[r+j].getOmega())  
            pose.set_omega((r+j+1), gen.getOmega())
            pose.set_phi((r+j+1), gen.getPhi())
            pose.set_psi((r+j+1), gen.getPsi())
            genes[r+j] = gen

        trial.setGenes(genes)         
        trial.setPose(pose)
        trial.setFitness(self.rosetta_conf.score3(pose))
        self.avaliacoes -= 1
        self.n_aval += 1
        return trial
        


    def gera_trial(self, seed, F_new, CR_new):
        
        trial = Individuo(F_new, CR_new)
        pose = self.rosetta_conf.get_generalPose()
        genes = []

        # loop para gerar novo indivíduo
        for j in range(0, self.D):
            r = random.uniform(0,1)
            gen = Gene()
                    
            if r <= trial.getCR() or j == random.randint(0, self.D):
                #phi, psi, omega = self.gera_mutacao(seed, F_new, seedA, indA, j)
                gen.setPhi(seed.genes[j].getPhi())
                gen.setPsi(seed.genes[j].getPsi())
                gen.setOmega(seed.genes[j].getOmega())
            else:
                gen.setPhi(self.target.genes[j].getPhi())
                gen.setPsi(self.target.genes[j].getPsi())
                gen.setOmega(self.target.genes[j].getOmega())
 
            # seta o objeto pose para calcular o fitness do individuo
            pose.set_omega((j+1), gen.getOmega())
            pose.set_phi((j+1), gen.getPhi())
            pose.set_psi((j+1), gen.getPsi())

            genes.append(gen)
                    
        trial.setGenes(genes)                
        trial.setPose(pose)
        trial.setFitness(self.rosetta_conf.score3(pose))
        self.avaliacoes -= 1
        self.n_aval += 1
        return trial


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


    # Retorna index de outra especie
    def get_outra_especie(self, tam, atual):
        r = random.randint(0, tam-1)
        while r == atual:
            r = random.randint(0, tam-1)
        return r

    # Retorna o index de um indivíduo da mesma espécie, diferente de i e 0
    def get_outro_ind(self, tam, i):
        r = random.randint(0, tam-1)
        while r == i or r == 0:
            r = random.randint(0, tam-1)
        return r


    # Retorna 1 indivíduo aleatório de toda a população, que não seja da especie atual
    def get_ind_qualquer(self, i):
        espec1 = random.randint(0, len(self.especies)-1)
        while espec1 == i:
            espec1 = random.randint(0, len(self.especies)-1)

        ind1 = random.randint(0, len(self.especies[espec1])-1)
        return self.especies[espec1][ind1]

    def forced_insertion_frag(self, ind):
        #score_before = ind.getFitness()
        #print('score antes %f \n' %score_before)
        pose = ind.getPose()

        self.rosetta_conf.get_9mer().apply(pose)
        ind = self.update_angle_from_pose(ind, pose)
        ind.setFitness(self.rosetta_conf.score3(pose))

        self.avaliacoes -= 1
        self.n_aval += 1

        #score_after = ind.getFitness()
        #print('score depois %f \n' %score_after)

        return ind


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
                print('%i %f' %(self.especies[i][j].getIndex(), self.especies[i][j].getFitness()))
            break

    def reinicia_populacao(self, ind):
        self.pop[:] = []
        for i in range(0, self.NP):
            novo = copy(ind)
            novo = self.forced_insertion_frag(novo)
            self.pop.append(novo)   

    def reinicia_especies(self, best):
        for i in range(0, len(self.especies)):
            tam = len(self.especies[i])
            corte = (int)(tam/3)
            while corte > 0:
                tam -= 1
                corte -=1
                ind = copy(best)
                index = self.especies[i][tam].getIndex()
               # ind = Individuo(self.F, self.CR)
               # ind = self.new_indiv_frag(ind)
                #ind = self.forced_insertion_frag(ind)
                ind = self.fragment_insert(ind=ind, n=100, mode=self.get_mode())
                ind.setIndex(index)
               # ind.setFitness(self.rosetta_conf.score3(ind.getPose()))
                self.especies[i][tam] = ind


    def otimiza(self):
        self.init_pop() # gera a população inicial   
        #self.printPop()
        #self.gera_vetor_dist()
        #self.printDistancias()
        #self.gera_especies()
        #self.printEspecies()  
        #print('especies %i \n'%len(self.especies))  
        print('populaçao %i \n' %len(self.pop)) 
        
        especAux = []
        m = 0
        n_ger = 1
        cont_best = 0
        index_seeds = []
                    
        while m < self.maxIteractions and self.avaliacoes > 0:  
            index_seeds[:] = []
            #if n_ger >= 5:
            n_ger = 0
             #   print('NOVA ESPECIE\n\n\n\n')
            self.gera_vetor_dist()
            self.gera_especies()     
            

            if cont_best >= 5:
                cont_best = 0
                self.reinicia_especies(copy(self.especies[0][0]))

            self.n_aval = 0       # numero de avaliações desde a geração da especie
            self.printEspecies()

            self.pop[:] = []
            best_index = 0    # armazena o índice do melhor indivíduo da população geral
            bad_index = 0     # armazena o índice do pior indivíduo da população geral
            soma = 0          # serve para calcular a média da população
            cont = 0          # conta o número de indivíduos na população geral

            
            for i in range(0, len(self.especies)):

                best_espec = 0                          # armazena o indice do best da especie
                seed = self.especies[i][0]              # semente da espécie atual
                aux = []                                # inicializa uma especie

                for j in range(0, len(self.especies[i])):

                    #if random.uniform(0,1) <= 0.02:
                    #    self.especies[i][j] = self.forced_insertion_frag(self.especies[i][j])

                    self.target = self.especies[i][j]   # pega o indivíduo atual

                    # Auto ajuste dos parâmetros jDE
                    F_new, CR_new = self.de_service.autoAjuste(self.target.getF(), self.target.getCR())

                    # Retorna outra semente
                    seedA = self.especies[self.get_outra_especie(len(self.especies), i)][0]

                    # Retorna individuo da mesma especie diferente de target e seed
                    indA = self.especies[i][self.get_outro_ind(len(self.especies[i]), j)]
                    ind1 = self.especies[i][self.get_outro_ind(len(self.especies[i]), j)]
                    ind2 = self.especies[i][self.get_outro_ind(len(self.especies[i]), j)]
                    ind3 = self.especies[i][self.get_outro_ind(len(self.especies[i]), j)]
                    ind4 = self.especies[i][self.get_outro_ind(len(self.especies[i]), j)]

                    # Retorna um indivíduo aleatório da população diferente da especie atual
                    indB = self.get_ind_qualquer(i)
                    indC = self.get_ind_qualquer(i)

                ###### CRUZAMENTO 

                    if self.mutation == 3:
                        if j != 0:
                            U = self.gera_trial_3(seed, F_new, CR_new)
                            U = self.fragment_insert(ind=copy(U), mode=self.get_mode())
                        else:
                            U = self.fragment_insert(ind=copy(seed), mode=self.get_mode())

                    elif self.mutation == 1:
                        U1 = self.gera_trial_2(seed, indA, ind1, F_new, CR_new)
                        if random.uniform(0,1) <= U1.getF():
                            U1 = self.fragment_insert(ind=copy(U1), mode=self.get_mode())
                        U2 = self.gera_trial_2(ind2, ind3, ind4, F_new, CR_new)
                        if random.uniform(0,1) <= U2.getF():
                            U2 = self.fragment_insert(ind=copy(U2), mode=self.get_mode())

                        if U1.getFitness() <= U2.getFitness():
                            U = U1
                        else:
                            U = U2

                    else:
                        # Gera vetor seed-to-seed
                        if j != 0:
                            U1 = self.gera_trial(seed, F_new, CR_new)
                            # MUTAÇÃO
                            if random.uniform(0,1) <= U1.getF():
                                U1 = self.fragment_insert(ind=copy(U1), mode=self.get_mode())
                        else:
                            U1 = self.gera_trial(seedA, F_new, CR_new)
                            # MUTAÇÃO
                            if random.uniform(0,1) <= U1.getF():
                                U1 = self.fragment_insert(ind=copy(U1), mode=self.get_mode())

                        # Gera vetor seed-to-rand
                        U2 = self.gera_trial(indB, F_new, CR_new)
                        # MUTAÇÃO
                        if random.uniform(0,1) <= U2.getF():
                            U2 = self.fragment_insert(ind=U2, mode=self.get_mode())

                        # COMPETIÇÃO ENTRE INTENSIFICAÇÃO E DIVERSIFICAÇÃO
                        if U1.getFitness() <= U2.getFitness():
                            U = U1                    
                        else:
                            U = U2

                ###### SELEÇÃO
                    if U.getFitness() <= self.target.getFitness():
                        self.pop.append(U)
                        aux.append(U)
                    else:
                        self.pop.append(self.target)
                        aux.append(self.target)
                    self.pop[cont].setIndex(cont)       # seta o index do indivíduo (posição na população)
                    aux[j].setIndex(cont)

                    # Verifica o melhor indivíduo de toda a população 
                    if self.pop[cont].getFitness() <= self.pop[best_index].getFitness():
                        best_index = cont
                    if self.pop[cont].getFitness() >= self.pop[bad_index].getFitness():
                        bad_index = cont

                    # Verifica o melhor indivíduo da espécie atual
                    if aux[j].getFitness() <= aux[best_espec].getFitness():
                        best_espec = j
                    
                    soma = soma + self.pop[cont].getFitness()
                    cont = cont+1       # incrementa o número de indivíduos

                # Atualiza a semente da espécie atual
                if best_espec != 0:
                    temp = aux[0]
                    aux[0] = aux[best_espec]
                    aux[best_espec] = temp
                especAux.append(aux)
                index_seeds.append(aux[0].getIndex())
            
            
            Log().best_score3(self.caminho+self.nome, self.pop[best_index].getFitness())
            Log().media_score3(self.caminho+self.nome, soma / self.NP)
            if self.pop[best_index].getFitness() >= self.min_fit:
                cont_best += 1
            else:
                cont_best = 0

            self.min_fit = self.pop[best_index].getFitness()      
            self.max_fit = self.pop[bad_index].getFitness()
                
            
            self.especies[:] = []
            self.especies = copy(especAux)
            especAux[:] = []
            m = m+1    # incrementa o número de gerações
            n_ger += 1

        for k in range(0, len(index_seeds)):
            rmsd = self.rosetta_conf.get_rmsd_from_native(self.pop[index_seeds[k]].getPose())
            print('index %i rmsd %f \n' %(index_seeds[k], rmsd))

        rmsd_best_pop = self.rosetta_conf.get_rmsd_from_native(self.pop[best_index].getPose())
        print('index best pop %i rmsd best pop %f \n' %(best_index,rmsd_best_pop))

        Log().best_score3(self.caminho+self.nome, self.pop[best_index].getFitness())
        Log().media_score3(self.caminho+self.nome, soma / self.NP)
        print('---------------------------------SCORE3 %f' %self.pop[best_index].getFitness())
        self.printPDB(self.pop[best_index].getPose())

        return self.pop[best_index].getPose()

                    

