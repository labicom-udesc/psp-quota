import os
import re


class ProteinLoader:
    def __init__(self, name=None, protocolo=None):
        self.name = name
        self.protocolo = protocolo

        self.fragset3_path = None
        self.fragset9_path = None

        self.fasta = None
        self.ssq8 = None      # SS prevista pelo mufold para Q8
        self.ssq3 = None      # SS prevista pelo mufold para Q3
        self.mufold_prob = [] # matriz com as probabilidades dos tipos SS
        self.psi_pred = None

        self.native_path = None

        self.protein_data_path = '../../protein_data'

    def load(self, name=None, protocolo=None):
        if name is None:
            name = self.name
        else:
            self.name = name

        if protocolo is None:
            protocolo = self.protocolo
        else:
            self.protocolo = protocolo

        if name is None and self.name is None:
            raise ValueError('No protein was specified')

        if protocolo is None and self.protocolo is None:
            raise ValueError('No protocol was specified')

        original = os.path.dirname(os.path.realpath(__file__))
        self.original = os.getcwd()

        base = self.protein_data_path

        dest = base + '/' + self.name

        if os.path.isdir(dest):
            os.chdir(dest)
            base_file = name + '.pdb'
            native_path = os.path.join(base, name, base_file)
            if os.path.isfile(native_path):
                self.native_path = native_path
            else:
                raise FileNotFoundError('Could not find %s' % native_path)

            if os.path.isfile(name + '.fasta'):
                with open(name + '.fasta', 'rt') as f:
                    self.fasta = ''
                    for line in f.readlines()[1:]:
                        l = line.rstrip()
                        self.fasta += l

                self.fasta = self.fasta.split('>')[0]

            if os.path.isfile(name + '.psipred.ss2'):
                with open(name + '.psipred.ss2', 'rt') as f:
                    self.psi_pred = ''
                    for line in f.readlines()[1:]:
                        tokens = re.sub(r"\s+", " ", line.rstrip().lstrip()).split(' ')
                        if len(tokens) < 3:
                            continue

                        self.psi_pred += tokens[2]


            if os.path.isfile('mufold/query.ss'):
                with open('mufold/query.ss', 'rt') as f:
                    self.ssq3 = ''
                    self.ssq8 = ''
                    seq = f.readlines()
                    self.ssq3 = seq[1].strip()
                    self.ssq8 = seq[3].strip()


            if os.path.isfile('mufold/query.fasta.q8prob'):
                with open('mufold/query.fasta.q8prob', 'rt') as f:
                    c = 0
                    for line in f.readlines():
                        if c == 0:
                            c += 1
                        else:
                            entries = re.split(',', line.strip())
                            self.mufold_prob.append(entries[0:7])


            f3 = base + '/' + name + '/output_' + protocolo + '/' + name + '.200.3mers'
            if os.path.isfile(f3):
                self.fragset3_path = f3
            else:
                raise FileNotFoundError('Could not find %s' % f3)

            f9 = base + '/' + name + '/output_' + protocolo + '/' + name + '.200.9mers'
            if os.path.isfile(f9):
                self.fragset9_path = f9
            else:
                raise FileNotFoundError('Could not find %s' % f9)
        else:
            raise FileNotFoundError('Could not find base folder: %s' % dest)

        os.chdir(original)

    def get_data(self):
        return (self.name, self.fasta, self.psi_pred, self.ssq3, self.ssq8, self.native_path, 
                self.fragset3_path, self.fragset9_path, self.mufold_prob)
