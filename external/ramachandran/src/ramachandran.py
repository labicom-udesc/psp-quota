import re
import os.path
import math
import random
import numpy as np
import os


class Ramachandran():
    def __init__(self):
        self.data = {}
        self.parsed_data = {}
        self.aa = {}

        self.aa['g'] = 'gly'
        self.aa['a'] = 'ala'
        self.aa['p'] = 'pro'
        self.aa['s'] = 'ser'
        self.aa['c'] = 'cys'
        self.aa['t'] = "thr"
        self.aa['v'] = 'val'
        self.aa['i'] = 'ile'
        self.aa['l'] = 'leu'
        self.aa['d'] = 'asp'
        self.aa['n'] = 'asn'
        self.aa['f'] = 'phe'
        self.aa['y'] = 'tyr'
        self.aa['h'] = 'his'
        self.aa['w'] = 'trp'
        self.aa['m'] = 'met'
        self.aa['e'] = 'glu'
        self.aa['q'] = 'gln'
        self.aa['k'] = 'lys'
        self.aa['r'] = 'arg'

    def load(self, path=""):
        if os.path.isfile(path):
            with open(path) as f:
                for line in f.readlines():
                    line = line.rstrip('\r\n')
                    tokens = re.sub("\s+", " ", line).split(' ')
                    if tokens[0][0] != "#":
                        name, phi, psi, prob, _ = tokens
                        phi = round(float((phi)))
                        psi = round(float((psi)))
                        prob = float(prob)
                        if name not in self.data:
                            a = {}
                            a[(phi, psi)] = prob
                            self.data[name] = a
                        else:
                            self.data[name][(phi, psi)] = prob
        else:
            if path == "":
                raise Exception("path is empty, please use keyword path to point to the file")
            else:
                raise Exception("Could not find: %s" % path)

    def process(self):
        f = None
        size = None
        dx = None
        for k, v in self.data.items():
            if f is None:
                f = len(v)
            else:
                if len(v) != f:
                    raise Exception("Uneven number of information")

        if round(math.sqrt(f)) == math.sqrt(f):
            size = round(math.sqrt(f))
            dx = (360) / size
        else:
            raise Exception("Unable to determine grid size")

        for k, v in self.data.items():
            self.parsed_data[k] = np.zeros((size, size))
            for a, b in self.data[k].items():
                phi, psi = a
                x = round((phi + 180) / dx)
                y = round((psi + 180) / dx)
                self.parsed_data[k][x, y] = self.data[k][a]

        self.size = size
        self.dx = dx

    def index(self, name, phi, psi):
        x = round((phi + 180) / self.dx)
        y = round((psi + 180) / self.dx)
        return self.parsed_data[name][x, y]

    def gnuprint(self, name, path=None):
        if path is None:
            for phi in range(-180, 180, 10):
                for psi in range(-180, 180, 10):
                    print("%4d %4d %.4f" % (phi, psi, self.index(name, phi, psi)))
                print()
        else:
            with open(path, "w") as f:
                for phi in range(-180, 180, 10):
                    for psi in range(-180, 180, 10):
                        f.write("%4d %4d %.4f\n" % (phi, psi, self.index(name, phi, psi)))
                    f.write("\n")

    def get_random_angle(self, name):
        if len(name) != 3:
            name = self.convert(name)

        r = random.random()
        c = 0.0
        x, y = 0, 0

        while c < r and x < self.size and y < self.size:
            x += 1
            if x == self.size:
                x = 0
                y += 1

            if y == 36:
                y -= 1

            c += self.parsed_data[name][x, y]

        phi = (x * self.dx) - 180 + random.random() * self.dx
        psi = (y * self.dx) - 180 + random.random() * self.dx

        return (phi, psi)

    def convert(self, name):
        return self.aa[name.lower()].lower()


if __name__ == "__main__":
    rama = Ramachandran()
    p1 = "/home/h3nnn4n/PyRosetta4.Release.python35.linux.release-147/setup/pyrosetta/database/scoring/score_functions/rama/shapovalov/kappa75/all.ramaProb"
    p2 = "/usr/local/lib/python3.5/dist-packages/pyrosetta-4.0-py3.5-linux-x86_64.egg/pyrosetta/database/scoring/score_functions/rama/shapovalov/kappa75/all.ramaProb"
    p3 = "/usr/local/lib/python3.6/site-packages/pyrosetta-2018.22+release.99b36feae43-py3.6-macosx-10.13-x86_64.egg/pyrosetta/database/scoring/score_functions/rama/shapovalov/kappa75/all.ramaProb"

    if os.path.exists(p1):
        rama.load(path=p1)
    elif os.path.exists(p2):
        rama.load(path=p2)
    else:
        rama.load(path=p3)

    rama.process()

    print(rama.get_random_angle('T'))

    import matplotlib.pyplot as plt
    for k, v in rama.parsed_data.items():
        plt.imshow(v.transpose(), cmap='viridis', interpolation='nearest', origin="lower")
        plt.savefig(k + ".png")
