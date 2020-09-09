import pyrosetta
import protein_loader
import os
import sys

sys.path.append('../external')

import ramachandran.src.ramachandran as rama  # pylint: disable=E0401

class RosettaConfig:

    def __init__(self, name):
        self.scorefxn = pyrosetta.get_fa_scorefxn()
        self.score0 = pyrosetta.create_score_function('score0')
        self.score3 = pyrosetta.create_score_function('score3')

        self.allatom = False
        
        self.allatom_switch = pyrosetta.SwitchResidueTypeSetMover('fa_standard')
        self.centroid_switch = pyrosetta.SwitchResidueTypeSetMover('centroid')

        self.fragset3 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(3)
        self.fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
        
        self.protein_loader = protein_loader.ProteinLoader()
        self.protein_loader.load(name)
        self.name, self.fasta, self.psi_pred, self.ssq3, self.ssq8, self.native_path, self.fragset3_path, self.fragset9_path, self.mufold_prob = self.protein_loader.get_data()

        self.native = pyrosetta.pose_from_pdb(self.native_path)

        self.dssp = pyrosetta.rosetta.protocols.moves.DsspMover()

        # ARRUMAR O PATH DO RAMACHANDRAN
        self.ramachandran = rama.Ramachandran()
        p1 = "/usr/local/lib/python3.6/dist-packages/pyrosetta-2019.39+release.93456a567a8-py3.6-linux-x86_64.egg/pyrosetta/database/scoring/score_functions/" + \
             "rama/shapovalov/kappa75/all.ramaProb"
        p2 = "/usr/lib/python3.5/site-packages/pyrosetta-4.0-py3.5-linux-x86_64.egg/pyrosetta/database/scoring/score_functions/rama/s" + \
             "hapovalov/kappa75/all.ramaProb"
        p3 = "/usr/local/lib/python3.6/dist-packages/pyrosetta-2019.39+release.93456a567a8-py3.6-linux-x86_64.egg/pyrosetta/da" + \
             "tabase/scoring/score_functions/rama/shapovalov/kappa75/all.ramaProb"

        if os.path.exists(p1):
            self.ramachandran.load(path=p1)
        elif os.path.exists(p2):
            self.ramachandran.load(path=p2)
        else:
            self.ramachandran.load(path=p3)

        self.ramachandran.process()

        self.movemap = pyrosetta.MoveMap()
        self.movemap.set_bb(True)

        cost = pyrosetta.rosetta.protocols.simple_moves.GunnCost()
        
        self.fragset3.read_fragment_file(self.fragset3_path)
        self.fragset9.read_fragment_file(self.fragset9_path)

        self.fast_relax = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefxn)

        self.mover_3mer = pyrosetta.rosetta.protocols.simple_moves.ClassicFragmentMover(self.fragset3, self.movemap)
        self.mover_9mer = pyrosetta.rosetta.protocols.simple_moves.ClassicFragmentMover(self.fragset9, self.movemap)
        self.mover_3mer_smooth = pyrosetta.rosetta.protocols.simple_moves.SmoothFragmentMover(self.fragset3, self.movemap, cost)
        self.mover_9mer_smooth = pyrosetta.rosetta.protocols.simple_moves.SmoothFragmentMover(self.fragset9, self.movemap, cost)

        self.mover_3mer.set_movemap(self.movemap)
        self.mover_9mer.set_movemap(self.movemap)
        self.mover_3mer_smooth.set_movemap(self.movemap)
        self.mover_9mer_smooth.set_movemap(self.movemap)

            
    def get_fast_relax(self):
        return self.fast_relax


    def copy_pose_to_allatom(self, pose):
        ap = pyrosetta.Pose()
        ap.assign(pose)

        if ap.is_centroid():
            self.allatom_switch.apply(ap)

        return ap


    def get_9mer(self):
        return self.mover_9mer

    def get_3mer(self):
        return self.mover_3mer

    def get_9mer_smooth(self):
        return self.mover_9mer_smooth

    def get_3mer_smooth(self):
        return self.mover_3mer_smooth

    def get_mer(self, mode):
        if mode == '3':
            return self.get_3mer()
        elif mode == '3s':
            return self.get_3mer_smooth()
        elif mode == '9':
            return self.get_9mer()
        elif mode == '9s':
            return self.get_9mer_smooth()
        else:
            return self.get_9mer()

    def get_new_mc(self, pose, score, temp=2.0):
        return pyrosetta.MonteCarlo(pose, score, temp)

    def get_new_seq_mover(self):
        return pyrosetta.rosetta.protocols.moves.SequenceMover()

    def get_new_trial_mover(self, seq, mc):
        return pyrosetta.TrialMover(seq, mc)

    def get_new_rep_mover(self, trial, n):
        return pyrosetta.rosetta.protocols.moves.RepeatMover(trial, n)

    def random_rama_angles_to_pose(self, pose):
        for k, v in enumerate(self.fasta):
            phi, psi = self.ramachandran.get_random_angle(v)
            pose.set_phi(k + 1, phi)
            pose.set_psi(k + 1, psi)

    def get_new_pose(self):
        pose = pyrosetta.pose_from_sequence(self.fasta)
        self.random_rama_angles_to_pose(pose)
        self.centroid_switch.apply(pose)
        return pose

    def get_generalPose(self):
        pose = pyrosetta.pose_from_sequence(self.fasta)
        self.centroid_switch.apply(pose)
        return pose

    def get_rmsd_from_pose(self, pose=None, pose2=None):
        if pose2 is None:
            raise TypeError('Rosetta pack with internal states is no longer suported')
        else:
            return pyrosetta.rosetta.core.scoring.CA_rmsd(pose, pose2)

    def get_rmsd_from_native(self, pose):
        return self.get_rmsd_from_pose(self.native, pose)