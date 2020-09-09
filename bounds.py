import random


class Bounds:
    def __init__(self):
        self.sdd = {}
        self.sdd["G"] = 0
        self.sdd["A"] = 0
        self.sdd["P"] = 0
        self.sdd["S"] = 1
        self.sdd["C"] = 1
        self.sdd["T"] = 1
        self.sdd["V"] = 1
        self.sdd["I"] = 2
        self.sdd["L"] = 2
        self.sdd["D"] = 2
        self.sdd["N"] = 2
        self.sdd["F"] = 2
        self.sdd["Y"] = 2
        self.sdd["H"] = 2
        self.sdd["W"] = 2
        self.sdd["M"] = 3
        self.sdd["E"] = 3
        self.sdd["Q"] = 3
        self.sdd["K"] = 4
        self.sdd["R"] = 4

        self.sdd["GLY"] = 0
        self.sdd["ALA"] = 0
        self.sdd["PRO"] = 0
        self.sdd["SER"] = 1
        self.sdd["CYS"] = 1
        self.sdd["THR"] = 1
        self.sdd["VAL"] = 1
        self.sdd["ILE"] = 2
        self.sdd["LEU"] = 2
        self.sdd["ASP"] = 2
        self.sdd["ASN"] = 2
        self.sdd["PHE"] = 2
        self.sdd["TYR"] = 2
        self.sdd["HIS"] = 2
        self.sdd["TRP"] = 2
        self.sdd["MET"] = 3
        self.sdd["GLU"] = 3
        self.sdd["GLN"] = 3
        self.sdd["LYS"] = 4
        self.sdd["ARG"] = 4

    def generateRandomAngles(self, ss):
        phi = random.uniform(self.getSecondaryLowerBound(ss, 0), self.getSecondaryUpperBound(ss, 0))
        psi = random.uniform(self.getSecondaryLowerBound(ss, 1), self.getSecondaryUpperBound(ss, 1))
        omega = 180.0 + random.gauss(0, 5.0)

        if omega > 180:
            omega -= 360

        return phi, psi, omega

    def generateRandomSidechainAngles(self, aa):
        r = []
        for i in range(self.getNumSideChainAngles(aa)):
            v = random.uniform(self.getSideChainLowerBound(i + 2, aa), self.getSideChainUpperBound(i + 2, aa))
            r.append(v)

        return r

    def getNumSideChainAngles(self, name):
        # if name not in self.sdd.keys():
        #     return 0
        # return self.sdd[name]
        if (name == 'GLY' or name == 'G') or (name == 'ALA' or name == 'A') or \
                (name == 'PRO' or name == 'P'):
                    return 0
        elif (name == 'SER' or name == 'S') or (name == 'CYS' or name == 'C') or \
                (name == "THR" or name == 'T') or (name == 'VAL' or name == 'V'):
                    return 1
        elif (name == 'ILE' or name == 'I') or (name == 'LEU' or name == 'L') or \
                (name == 'ASP' or name == 'D') or (name == 'ASN' or name == 'N') or \
                (name == 'PHE' or name == 'F') or (name == 'TYR' or name == 'Y') or \
                (name == 'HIS' or name == 'H') or (name == 'TRP' or name == 'W'):
                    return 2
        elif (name == 'MET' or name == 'M') or (name == 'GLU' or name == 'E') or \
                (name == 'GLN' or name == 'Q'):
                    return 3
        elif (name == 'LYS' or name == 'K') or (name == 'ARG' or name == 'R'):
            return 4
        return 0

    def getSecondaryLowerBound(self, secondaryStructure, i):
        if secondaryStructure == "H":
            if i == 0:
                return -67
            elif i == 1:
                return -57
        elif secondaryStructure == "B" or secondaryStructure == "E":
            if i == 0:
                return -130
            elif i == 1:
                return 110
        elif secondaryStructure == "G":
            if i == 0:
                return -59
            elif i == 1:
                return -36
        elif secondaryStructure == "I":
            if i == 0:
                return -67
            elif i == 1:
                return -80
        elif secondaryStructure == "T" or secondaryStructure == "S" or secondaryStructure == "C":
            return -180
        else:
            raise Exception("ERROR! Recondary Chain not Found! ", secondaryStructure, i)

    def getSecondaryUpperBound(self, secondaryStructure, i):
        if secondaryStructure == "H":
            if i == 0:
                return -47
            elif i == 1:
                return -37
        elif secondaryStructure == "B" or secondaryStructure == "E":
            if i == 0:
                return -110
            elif i == 1:
                return 130
        elif secondaryStructure == "G":
            if i == 0:
                return -39
            elif i == 1:
                return 16
        elif secondaryStructure == "I":
            if i == 0:
                return -47
            elif i == 1:
                return -60
        elif secondaryStructure == "T" or secondaryStructure == "S" or secondaryStructure == "C":
            return 180
        else:
            raise Exception("ERROR! Recondary Chain not Found! ", secondaryStructure, i)

    def getSideChainLowerBound(self, i, name):
        if name == "ARG" or name == "R":
            if i == 2:
                return -177
            elif i == 3:
                return -167
            elif i == 4:
                return -65
            elif i == 5:
                return -175
        elif name == "LYS" or name == "K":
            if i == 2:
                return -177
            elif i == 3:
                return -68
            elif i == 4:
                return -68
            elif i == 5:
                return -65
        elif name == "MET" or name == "M":
            if i == 2:
                return -177
            elif i == 3:
                return -65
            elif i == 4:
                return -75
        elif name == "GLU" or name == "E":
            if i == 2:
                return -177
            elif i == 3:
                return -80
            elif i == 4:
                return -60
        elif name == "GLN" or name == "Q":
            if i == 2:
                return -177
            elif i == 3:
                return -75
            elif i == 4:
                return -100
        elif name == "ASP" or name == "D":
            if i == 2:
                return -177
            elif i == 3:
                return -60
        elif name == "ASN" or name == "N":
            if i == 2:
                return -177
            elif i == 3:
                return -80
        elif name == "ILE" or name == "I":
            if i == 2:
                return -177
            elif i == 3:
                return -60
        elif name == "LEU" or name == "L":
            if i == 2:
                return -177
            elif i == 3:
                return 65
        elif name == "HIS" or name == "H":
            if i == 2:
                return -177
            elif i == 3:
                return -165
        elif name == "TRP" or name == "W":
            if i == 2:
                return -177
            elif i == 3:
                return -105
        elif name == "TYR" or name == "Y":
            if i == 2:
                return -177
            elif i == 3:
                return -85
        elif name == "PHE" or name == "F":
            if i == 2:
                return -177
            elif i == 3:
                return -85
        elif name == "THR" or name == "T":
            if i == 2:
                return -177
        elif name == "VAL" or name == "V":
            if i == 2:
                return -60
        elif name == "SER" or name == "S":
            if i == 2:
                return -177
        elif name == "CYS" or name == "C":
            if i == 2:
                return -177

        print('PANIC, residue not found!!', 1, name)

    def getSideChainUpperBound(self, i, name):
        if name == "ARG" or name == "R":
            if i == 2:
                return 62
            elif i == 3:
                return 180
            elif i == 4:
                return 180
            elif i == 5:
                return 180
        elif name == "LYS" or name == "K":
            if i == 2:
                return 62
            elif i == 3:
                return 180
            elif i == 4:
                return 180
            elif i == 5:
                return 180
        elif name == "MET" or name == "M":
            if i == 2:
                return 62
            elif i == 3:
                return 180
            elif i == 4:
                return 180
        elif name == "GLU" or name == "E":
            if i == 2:
                return 70
            elif i == 3:
                return 180
            elif i == 4:
                return 60
        elif name == "GLN" or name == "Q":
            if i == 2:
                return 70
            elif i == 3:
                return 180
            elif i == 4:
                return 100
        elif name == "ASP" or name == "D":
            if i == 2:
                return 62
            elif i == 3:
                return 65
        elif name == "ASN" or name == "N":
            if i == 2:
                return 62
            elif i == 3:
                return 120
        elif name == "ILE" or name == "I":
            if i == 2:
                return 62
            elif i == 3:
                return 170
        elif name == "LEU" or name == "L":
            if i == 2:
                return 62
            elif i == 3:
                return 175
        elif name == "HIS" or name == "H":
            if i == 2:
                return 62.0
            elif i == 3:
                return 165.0
        elif name == "TRP" or name == "W":
            if i == 2:
                return 62
            elif i == 3:
                return 95
        elif name == "TYR" or name == "Y":
            if i == 2:
                return 62
            elif i == 3:
                return 90
        elif name == "PHE" or name == "F":
            if i == 2:
                return 62
            elif i == 3:
                return 90
        elif name == "THR" or name == "T":
            if i == 2:
                return 62
        elif name == "VAL" or name == "V":
            if i == 2:
                return 175
        elif name == "SER" or name == "S":
            if i == 2:
                return 62
        elif name == "CYS" or name == "C":
            if i == 2:
                return 62

        print('PANIC, residue not found!!', 1, name)
