import pyrosetta
from rosetta_config import RosettaConfig
from jDE import jDE
from dsmDE import dsmDE
from dsm import dsm
from dsmDE2 import dsm2
from dsm_paralelo import dsm_paralelo
import yaml
from repacker import Repacker
from copy import copy



def main():

    pyrosetta.init()

    pose = pyrosetta.Pose()

    with open("config.yaml", 'r') as stream:
        try:
            config = yaml.load(stream)
            proteinName = config['proteinName']
        except yaml.YAMLError as exc:                
            print(exc)
    rosetta_conf = RosettaConfig(proteinName)
    
    #de = jDE(rosetta_conf)
    #dsm = dsmDE(rosetta_conf)
    #de = dsm(rosetta_conf)
    de = dsm2(rosetta_conf)
    #de = dsm_paralelo(rosetta_conf)
    pose.assign(de.otimiza())


    repack = Repacker(rosetta_conf)
    pose = repack.repack(pose)
    
    score  = rosetta_conf.scorefxn(pose)

    print('---------------------------------SCOREFXN %f' %score)

    pose.dump_pdb("teste_"+proteinName+"_repacked.pdb")
 
    

if __name__ == '__main__':
    main()
