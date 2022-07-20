from ..config_win import *
from utils import seed_all
import torch


def inference():
    print('\n')
    arg_helper = EquiBindArgs()
    args = arg_helper.get_all_args()
    seed_all(args.seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() and args.device == 'cuda' else "cpu")
    checkpoint = torch.load(args.checkpoint, map_location=device)
    model = None
    all_ligs_coords_corrected = []
    all_intersection_losses = []
    all_intersection_losses_untuned = []
    all_ligs_coords_pred_untuned = []
    all_ligs_coords = []
    all_ligs_keypts = []
    all_recs_keypts = []
    all_names = []
    dp = args.dataset_params
    use_rdkit_coords = args.use_rdkit_coords if args.use_rdkit_coords != None else args.dataset_params[
        'use_rdkit_coords']
    for name in protein_names:
        print('Processing {}:' + name)


def test():
    pass


if __name__ == '__main__':
    inference()
