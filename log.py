class Log:


    def best_score3(self, nome, valor):
        arquivo = open(nome+'_best.txt', 'a')
        arquivo.write('%f\n' %valor)
        arquivo.close()

    def rmsd_diversity(self, nome, rmsd, diversidade, psc):
        arquivo = open(nome+'_diversity.txt', 'a')
        arquivo.write('%f %f %f\n' %(rmsd, diversidade, psc))
        arquivo.close()

    def media_score3(self, nome, valor):
        arquivo = open(nome+'_media.txt', 'a')
        arquivo.write('%f\n' %valor)
        arquivo.close()

    def rmsd(self, nome, especie, rmsd, energia):
        arquivo = open(nome+'_rmsd.txt', 'a')
        arquivo.write('%i %f %f\n' %(especie, rmsd, energia))
        arquivo.close()

    def repack(self, caminho, nome, rmsd, gdt, energia, tempo):
        arquivo = open(caminho+'_repack.txt', 'a')
        arquivo.write('%s %f %f %f %f\n' %(nome, rmsd, gdt, energia, tempo))
        arquivo.close()

    def time(self, nome, tempo, valor):
        arquivo = open(nome+'_time.txt', 'a')
        arquivo.write('%s: %f\n' %(tempo, valor))
        arquivo.close()

    def mapa(self, nome, trial_fit, x1, trial_ss, target_fit, x2, target_ss, dist_trial, dist_target, v_Trial, v_Target):
        arquivo = open(nome+'_mapa.txt', 'a')
        arquivo.write('Trial fitness: %f\n' %(trial_fit))
        arquivo.write('Trial mapa: %f\n' %dist_trial)
        arquivo.write('Trial ss: %s\n' %trial_ss)
        arquivo.write('Trial ss: %i\n' %x1)
        
        arquivo.write('Target fitness: %f\n' %(target_fit))
        arquivo.write('Target mapa: %f\n' %dist_target)
        arquivo.write('Target ss: %s\n' %target_ss)
        arquivo.write('Target ss: %i\n' %x2)

        arquivo.write('Trial: ')
        for i in range(0, len(v_Trial)):
            arquivo.write('%f ' %v_Trial[i])
        arquivo.write('\n')

        arquivo.write('Target: ')
        for j in range(0, len(v_Target)):
            arquivo.write('%f ' %(v_Target[j]))
        arquivo.write('\n\n\n')

        arquivo.close()
  