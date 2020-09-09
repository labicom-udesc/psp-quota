class Repacker:

    def __init__(self, rosetta_conf):        
        self.rosetta_conf = rosetta_conf

    # Repack deve receber um parâmetro Pose para calcular os ângulos da cadeia lateral
    def repack(self, pose):
        repack = self.rosetta_conf.get_fast_relax()
        best = self.rosetta_conf.copy_pose_to_allatom(pose)
        repack.apply(best)        

        return best