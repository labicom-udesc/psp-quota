class Log:


    def best_score3(self, nome, aval, valor):
        arquivo = open(nome+'_best.txt', 'a')
        arquivo.write('%i %f\n' %(aval,valor))
        arquivo.close()


    def media_score3(self, nome, aval, valor):
        arquivo = open(nome+'_media.txt', 'a')
        arquivo.write('%i %f\n' %(aval, valor))
        arquivo.close()

    def rmsd(self, nome, especie, rmsd, energia):
        arquivo = open(nome+'_rmsd.txt', 'a')
        arquivo.write('%i %f %f\n' %(especie, rmsd, energia))
        arquivo.close()

    def repack(self, nome, rmsd, energia):
        arquivo = open(nome+'_repack.txt', 'a')
        arquivo.write('%f %f\n' %(rmsd, energia))
        arquivo.close()

    def time(self, nome, tempo, valor):
        arquivo = open(nome+'_time.txt', 'a')
        arquivo.write('%s: %f\n' %(tempo, valor))
        arquivo.close()
  