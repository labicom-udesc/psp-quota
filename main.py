import pyrosetta
from rosetta_config import RosettaConfig
#from jDE import jDE
#from dsmDE import dsmDE
#from dsm import dsm
from dsmDE2 import dsm2
#from dsm_paralelo import dsm_paralelo
import yaml
from repacker import Repacker
from copy import copy
from log import Log
import time
import random

def main():

    pyrosetta.init()

    pose = pyrosetta.Pose()

    with open("config.yaml", 'r') as stream:
        try:
            config = yaml.load(stream)
            proteinName = config['proteinName']
            protocolo = config['protocolo']
        except yaml.YAMLError as exc:                
            print(exc)
    rosetta_conf = RosettaConfig(proteinName, protocolo)
    caminho = 'testes_protocol/'+protocolo+'/'+proteinName+'/'
    nome = str(random.random())
    
    #de = jDE(rosetta_conf)
    #dsm = dsmDE(rosetta_conf)
    #de = dsm(rosetta_conf)
    time_repack_start = time.time()
    de = dsm2(rosetta_conf, caminho, nome)
    #de = dsm_paralelo(rosetta_conf)
    pose.assign(de.otimiza())

    repack = Repacker(rosetta_conf)
    pose = repack.repack(pose)
    
    score  = rosetta_conf.scorefxn(pose)

    time_repack_end = time.time()

    rmsd = rosetta_conf.get_rmsd_from_native(pose) 

    Log().repack(caminho+nome, rmsd, score)
    Log().time(caminho+nome, 'Tempo com repack', time_repack_end - time_repack_start)

    #print('---------------------------------SCOREFXN %f' %score)

    pose.dump_pdb(caminho+nome+"_repacked.pdb")
 
    

if __name__ == '__main__':
    main()
