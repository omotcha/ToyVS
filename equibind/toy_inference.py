import os.path

from configs.config_win import *
from utils import seed_all
import torch
from rdkit import Chem
from dbutil.dbutil import *
from process_mols import read_molecule, get_lig_graph_revised, \
    get_rec_graph, get_geometry_graph, get_geometry_graph_ring, \
    get_receptor_inference
from geometry_utils import *
from train import load_model
from losses import *
from copy import deepcopy
from rdkit.Geometry import Point3D
from rdkit.Chem import AllChem
import pandas as pd


def inference_and_score():
    print('\n')
    arg_helper = EquiBindArgs()
    args = arg_helper.get_all_args()
    seed_all(args['seed'])
    device = torch.device("cuda:0" if torch.cuda.is_available() and args['device'] == 'cuda' else "cpu")
    checkpoint = torch.load(args['checkpoint'], map_location=device)
    model = None
    all_ligs_coords_corrected = []
    all_intersection_losses = []
    all_intersection_losses_untuned = []
    all_ligs_coords_pred_untuned = []
    all_ligs_coords = []
    all_ligs_keypts = []
    all_recs_keypts = []
    all_names = []
    dp = args['dataset_params']
    use_rdkit_coords = args['use_rdkit_coords'] if args['use_rdkit_coords'] is not None \
        else args['dataset_params']['use_rdkit_coords']
    db_helper = DBUtil()
    num_ligs = db_helper.get_num_rows('smiles2k')
    num_ligs = 1
    for name in protein_names:
        print('Processing {}:'.format(name))
        rec_path = os.path.join(protein_path, '{}.pdb'.format(name))
        for j in range(num_ligs):
            lig = AllChem.AddHs(Chem.MolFromSmiles(db_helper.fetch_smiles_by_index(j, 'smiles2k')))
            AllChem.EmbedMolecule(lig, useExpTorsionAnglePrefs=False, useBasicKnowledge=False)
            rec, rec_coords, c_alpha_coords, n_coords, c_coords = get_receptor_inference(rec_path)
            rec_graph = get_rec_graph(rec, rec_coords, c_alpha_coords, n_coords, c_coords,
                                      use_rec_atoms=dp['use_rec_atoms'], rec_radius=dp['rec_graph_radius'],
                                      surface_max_neighbors=dp['surface_max_neighbors'],
                                      surface_graph_cutoff=dp['surface_graph_cutoff'],
                                      surface_mesh_cutoff=dp['surface_mesh_cutoff'],
                                      c_alpha_max_neighbors=dp['c_alpha_max_neighbors'])

            lig_graph = get_lig_graph_revised(lig, name, max_neighbors=dp['lig_max_neighbors'],
                                              use_rdkit_coords=use_rdkit_coords, radius=dp['lig_graph_radius'])
            if 'geometry_regularization' in dp and dp['geometry_regularization']:
                geometry_graph = get_geometry_graph(lig)
            elif 'geometry_regularization_ring' in dp and dp['geometry_regularization_ring']:
                geometry_graph = get_geometry_graph_ring(lig)
            else:
                geometry_graph = None
            start_lig_coords = lig_graph.ndata['x']
            # Randomly rotate and translate the ligand.
            rot_T, rot_b = random_rotation_translation(translation_distance=5)
            if use_rdkit_coords:
                lig_coords_to_move = lig_graph.ndata['new_x']
            else:
                lig_coords_to_move = lig_graph.ndata['x']
            mean_to_remove = lig_coords_to_move.mean(dim=0, keepdims=True)
            input_coords = (rot_T @ (lig_coords_to_move - mean_to_remove).T).T + rot_b
            lig_graph.ndata['new_x'] = input_coords

            if model is None:
                model = load_model(args, data_sample=(lig_graph, rec_graph), device=device)
                model.load_state_dict(checkpoint['model_state_dict'])
                model.to(device)
                model.eval()

            with torch.no_grad():
                geometry_graph = geometry_graph.to(device) if geometry_graph != None else None
                ligs_coords_pred_untuned, ligs_keypts, recs_keypts, rotations, translations, geom_eg_loss = model(
                    lig_graph.to(device), rec_graph.to(device), geometry_graph, complex_names=[name], epoch=0)

                for lig_coords_pred_untuned, lig_coords, lig_keypts, rec_keypts, rotation, translation in zip(
                        ligs_coords_pred_untuned, [start_lig_coords], ligs_keypts, recs_keypts, rotations,
                        translations, ):
                    all_intersection_losses_untuned.append(
                        compute_revised_intersection_loss(lig_coords_pred_untuned.detach().cpu(), rec_graph.ndata['x'],
                                                          alpha=0.2, beta=8, aggression=0))
                    all_ligs_coords_pred_untuned.append(lig_coords_pred_untuned.detach().cpu())
                    all_ligs_coords.append(lig_coords.detach().cpu())
                    all_ligs_keypts.append(((rotation @ (lig_keypts).T).T + translation).detach().cpu())
                    all_recs_keypts.append(rec_keypts.detach().cpu())
                if args['run_corrections']:
                    prediction = ligs_coords_pred_untuned[0].detach().cpu()
                    lig_input = deepcopy(lig)
                    conf = lig_input.GetConformer()
                    for i in range(lig_input.GetNumAtoms()):
                        x, y, z = input_coords.numpy()[i]
                        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

                    lig_equibind = deepcopy(lig)
                    conf = lig_equibind.GetConformer()
                    for i in range(lig_equibind.GetNumAtoms()):
                        x, y, z = prediction.numpy()[i]
                        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

                    coords_pred = lig_equibind.GetConformer().GetPositions()

                    Z_pt_cloud = coords_pred
                    rotable_bonds = get_torsions([lig_input])
                    new_dihedrals = np.zeros(len(rotable_bonds))
                    for idx, r in enumerate(rotable_bonds):
                        new_dihedrals[idx] = get_dihedral_vonMises(lig_input, lig_input.GetConformer(), r, Z_pt_cloud)
                    optimized_mol = apply_changes(lig_input, new_dihedrals, rotable_bonds)

                    coords_pred_optimized = optimized_mol.GetConformer().GetPositions()
                    R, t = rigid_transform_Kabsch_3D(coords_pred_optimized.T, coords_pred.T)
                    coords_pred_optimized = (R @ (coords_pred_optimized).T).T + t.squeeze()
                    all_ligs_coords_corrected.append(coords_pred_optimized)

                    conf = optimized_mol.GetConformer()
                    for i in range(optimized_mol.GetNumAtoms()):
                        x, y, z = coords_pred_optimized[i]
                        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

                    # ecif::autogluon scoring
                    from ecif.util.ECIF import ECIF
                    ecif_helper = ECIF(2016)
                    ecif = ecif_helper.get_ecif(rec_path, optimized_mol, float(6.0))
                    ld = ecif_helper.get_ligand_features_by_mol(optimized_mol)
                    m = pickle.load(open(ecif_model, 'rb'))
                    data = ecif + list(ld)
                    cols = ecif_helper.get_possible_ecif() + ecif_helper.get_ligand_descriptors()
                    data_f = pd.DataFrame([data], columns=cols)
                    pred = m.predict(data_f)[0]
                    print("prediction: {}".format(pred))

                    # block_optimized = Chem.MolToMolBlock(optimized_mol)
                    # print('Writing predictions: ')
                    # with open(os.path.join(data_dir, 'equi_output', 'test.sdf'), "w") as newfile:
                    #     newfile.write(block_optimized)

    db_helper.__del__()


def test():
    pass


if __name__ == '__main__':
    inference_and_score()
