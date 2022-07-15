import pyrosetta
from rosetta_config import RosettaConfig
from dsmDE2 import dsm2
import yaml
from repacker import Repacker
from copy import copy
from log import Log
import time
import random

def main():

    pyrosetta.init("-mute all")

    pose = pyrosetta.Pose()

    with open("config.yaml", 'r') as stream:
        try:
            config = yaml.load(stream)
            proteinName = config['proteinName']
            protocolo = config['protocolo']
        except yaml.YAMLError as exc:                
            print(exc)
    rosetta_conf = RosettaConfig(proteinName, protocolo)
    caminho = 'testes/'+protocolo+'/'+proteinName+'/'
    nome = str(random.random())
    
    de = dsm2(rosetta_conf, caminho, nome)
    de.otimiza() 
    

if __name__ == '__main__':
    main()
