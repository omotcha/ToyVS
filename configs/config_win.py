"""
platform: win
env: any
name: config_win.py
configurations of ToyVS
"""

import os
import yaml

project_dir = "D:\\AlexYu\\ToyVS"
data_dir = os.path.join(project_dir, "data")
# in this configuration, 2 csv files with nearly 2k and 8m lines of smiles strings are used as data source
data_2k = os.path.join(data_dir, "smiles_2k.csv")
data_8m = os.path.join(data_dir, "smiles_8m.csv")
output_dir = os.path.join(data_dir, "equi_output")
log_dir = os.path.join(project_dir, 'logs')
protein_names = ['7AMA']
protein_path = os.path.join(data_dir, 'proteins')
num_process = 4
ecif_model = os.path.join(project_dir, 'ecif', 'models', 'ecif_gbt_6.0.pkl')
ecif_pdb_atom_keys = os.path.join(project_dir, 'ecif', 'ECIF_PDB_Atom_Keys.csv')


class EquiBindArgs:
    def __init__(self):
        self.device = 'cpu'
        self.run_corrections = True
        """
        generates the coordinates of the ligand with rdkit instead of using the provided conformer. 
        If you already have a 3D structure that you want to use as initial conformer, then leave this as False
        """
        self.use_rdkit_coords = False
        self.save_trajectories = False
        self.num_confs = False
        self.run_dir = os.path.join(project_dir, 'run', 'flexible_self_docking')
        self.checkpoint = os.path.join(self.run_dir, 'best_checkpoint.pt')
        self.metrics = ['pearsonr', 'rsquared', 'mean_rmsd', 'median_rmsd',
                        'median_centroid_distance', 'centroid_distance_less_than_2',
                        'mean_centroid_distance', 'kabsch_rmsd',
                        'rmsd_less_than_2', 'rmsd_less_than_5']
        self.other_args = {}
        self.arg_list = ['device', 'run_corrections', 'use_rdkit_coords',
                         'save_trajectories', 'num_confs', 'run_dir',
                         'checkpoint', 'metrics', 'other_args', 'arg_list']
        with open(os.path.join(os.path.dirname(self.checkpoint), 'train_arguments.yaml'), 'r') as arg_file:
            checkpoint_dict = yaml.load(arg_file, Loader=yaml.FullLoader)
        for key, value in checkpoint_dict.items():
            if key not in self.arg_list:
                if isinstance(value, list):
                    for v in value:
                        self.other_args[key].append(v)
                else:
                    self.other_args[key] = value

    def get_all_args(self):
        """
        get all the arguments needed for EquiBind input
        :return:
        """
        ret = self.other_args
        for arg in self.arg_list:
            ret[arg] = self.__getattribute__(arg)
        return ret


if __name__ == '__main__':
    helper = EquiBindArgs()
    all_args = helper.get_all_args()
    print(all_args)
