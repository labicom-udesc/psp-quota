#!/usr/bin/python3.5
# -*- coding: utf-8 -*-

import re
import os
import random
import os.path
import numpy as np


class APL():
    def __init__(self):
        self.data = {}
        self.ss = set()
        self.aa = set()

        self.size = 360
        self.dx = 360 / self.size

        self.names = {}

        self.names['g'] = 'GLY'
        self.names['a'] = 'ALA'
        self.names['p'] = 'PRO'
        self.names['s'] = 'SER'
        self.names['c'] = 'CYS'
        self.names['t'] = "THR"
        self.names['v'] = 'VAL'
        self.names['i'] = 'ILE'
        self.names['l'] = 'LEU'
        self.names['d'] = 'ASP'
        self.names['n'] = 'ASN'
        self.names['f'] = 'PHE'
        self.names['y'] = 'TYR'
        self.names['h'] = 'HIS'
        self.names['w'] = 'TRP'
        self.names['m'] = 'MET'
        self.names['e'] = 'GLU'
        self.names['q'] = 'GLN'
        self.names['k'] = 'LYS'
        self.names['r'] = 'ARG'

    def load(self, path="./apl"):
        dx = self.dx
        if os.path.isdir(path):
            for _, _, fnames in os.walk(path):
                for f in fnames:
                    sp = f.split('_')
                    self.aa.add(sp[0])
                    self.ss.add(sp[1])
                    key = (sp[0], sp[1])
                    self.data[key] = np.zeros((self.size, self.size))

                    p = path + "/" + f

                    if ".swp" not in p:
                        with open(p) as f:
                            for line in f.readlines():
                                line = line.rstrip('\r\n')
                                tokens = re.sub("\s+", " ", line).split(' ')

                                phi = float(tokens[0])
                                psi = float(tokens[1])
                                p = float(tokens[2])

                                x = round((phi + 180) / dx)
                                y = round((psi + 180) / dx)

                                self.data[key][x, y] = p

        else:
            if path == "":
                raise Exception("path is empty, please use keyword path to point to the file")
            else:
                raise Exception("Could not find: %s" % path)

    def get_map(self, aa, ss):
        if aa not in self.aa:
            raise Exception("> %s < not found in the Aminoacid list" % aa)

        if ss not in self.ss:
            raise Exception("> %s < not found in the secondary Structure list" % ss)

        k = (aa, ss)

        return self.data[k]

    def get_random_angle(self, aa, ss):
        if len(aa) != 3:
            aa = self.convert(aa)

        k = (aa, ss)

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

            c += self.data[k][x, y]

        phi = (x * self.dx) - 180 + random.random() * self.dx
        psi = (y * self.dx) - 180 + random.random() * self.dx

        return (phi, psi)

    def convert(self, name):
        return self.aa[name.lower()].upper()


if __name__ == "__main__":
    apl = APL()
    apl.load()

    # for _ in range(10):
    #     print(apl.get_random_angle('GLU', 'T'))

    import matplotlib.pyplot as plt
    for k, v in apl.data.items():
        print(k)
        plt.imshow(v.transpose(), cmap='viridis', interpolation='nearest', origin="lower")
        plt.savefig("%s_%s.png" % (k))
