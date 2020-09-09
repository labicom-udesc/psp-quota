class Log:


    def best_score3(self, nome, valor):
        arquivo = open(nome+'_best.txt', 'a')
        arquivo.write('%f\n' %valor)
        arquivo.close()


    def media_score3(self, nome, valor):
        arquivo = open(nome+'_media.txt', 'a')
        arquivo.write('%f\n' %valor)
        arquivo.close()
  