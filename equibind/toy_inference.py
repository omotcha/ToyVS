import multiprocessing
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
import time
import warnings
import multiprocessing as mp

warnings.filterwarnings("ignore")


def inference_and_score_mp(lig_id_list):
    """
    the multi-process worker for Equibind inferring and ECIF scoring
    :param lig_id_list: list of ligand id to be processed
    :return:
    """
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
    table_name = "distinct_smiles2k"
    ######################
    # only for test
    # table_name = "smiles2k"
    ######################
    m = pickle.load(open(ecif_model, 'rb'))
    for name in protein_names:
        # print('Processing {}:'.format(name))
        err_ids = []
        rec_path = os.path.join(protein_path, '{}.pdb'.format(name))
        rec, rec_coords, c_alpha_coords, n_coords, c_coords = get_receptor_inference(rec_path)
        rec_graph = get_rec_graph(rec, rec_coords, c_alpha_coords, n_coords, c_coords,
                                  use_rec_atoms=dp['use_rec_atoms'], rec_radius=dp['rec_graph_radius'],
                                  surface_max_neighbors=dp['surface_max_neighbors'],
                                  surface_graph_cutoff=dp['surface_graph_cutoff'],
                                  surface_mesh_cutoff=dp['surface_mesh_cutoff'],
                                  c_alpha_max_neighbors=dp['c_alpha_max_neighbors'])
        for j in lig_id_list:
            # for j in [89]:
            lig = AllChem.AddHs(
                Chem.MolFromSmiles(db_helper.fetch_canonical_smiles_by_index(j, table_name)))
            AllChem.EmbedMolecule(lig, useExpTorsionAnglePrefs=False, useBasicKnowledge=False)
            try:
                lig_graph = get_lig_graph_revised(
                    lig, name, max_neighbors=dp['lig_max_neighbors'],
                    use_rdkit_coords=use_rdkit_coords, radius=dp['lig_graph_radius'])
            except:
                err_ids.append(j)
                continue
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
                geometry_graph = geometry_graph.to(device) if geometry_graph is not None else None
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
                    # ecif::gbt scoring

                    from ecif.util.ECIF import ECIF
                    ecif_helper = ECIF(2016)
                    ecif_helper.cache_target(rec_path)
                    ecif = ecif_helper.get_ecif_cached(optimized_mol, float(6.0))
                    ld = ecif_helper.get_ligand_features_by_mol(optimized_mol)

                    data = ecif + list(ld)
                    cols = ecif_helper.get_possible_ecif() + ecif_helper.get_ligand_descriptors()
                    data_f = pd.DataFrame([data], columns=cols)
                    pred = m.predict(data_f)[0]
                    db_helper.insert_prediction(j, optimized_mol, pred, 'results2k')
        for err_id in err_ids:
            db_helper.insert_error(err_id, table_name)
    db_helper.__del__()
    return


def inference_and_score(table_name, lig_id_list):
    """
    Equibind inferring and ECIF scoring without multi-process
    :param table_name: table name
    :param lig_id_list: list of ligand id to be processed
    :return:
    """
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
    ######################
    # only for test
    # table_name = "smiles2k"
    ######################
    equibind_time = 0.0
    ecif_time = 0.0
    m = pickle.load(open(ecif_model, 'rb'))
    for name in protein_names:
        print('Processing {}:'.format(name))
        err_ids = []
        rec_path = os.path.join(protein_path, '{}.pdb'.format(name))
        rec, rec_coords, c_alpha_coords, n_coords, c_coords = get_receptor_inference(rec_path)
        rec_graph = get_rec_graph(rec, rec_coords, c_alpha_coords, n_coords, c_coords,
                                  use_rec_atoms=dp['use_rec_atoms'], rec_radius=dp['rec_graph_radius'],
                                  surface_max_neighbors=dp['surface_max_neighbors'],
                                  surface_graph_cutoff=dp['surface_graph_cutoff'],
                                  surface_mesh_cutoff=dp['surface_mesh_cutoff'],
                                  c_alpha_max_neighbors=dp['c_alpha_max_neighbors'])
        for j in lig_id_list:
            lig = AllChem.AddHs(
                Chem.MolFromSmiles(db_helper.fetch_canonical_smiles_by_index(j, table_name)))
            AllChem.EmbedMolecule(lig, useExpTorsionAnglePrefs=False, useBasicKnowledge=False)
            equibind_start = time.perf_counter()
            try:
                lig_graph = get_lig_graph_revised(
                    lig, name, max_neighbors=dp['lig_max_neighbors'],
                    use_rdkit_coords=use_rdkit_coords, radius=dp['lig_graph_radius'])
            except:
                err_ids.append(j)
                continue
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
                geometry_graph = geometry_graph.to(device) if geometry_graph is not None else None
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

                    equibind_end = time.perf_counter()
                    equibind_time = equibind_time + (equibind_end - equibind_start)
                    # ecif::gbt scoring

                    ecif_start = time.perf_counter()
                    from ecif.util.ECIF import ECIF
                    ecif_helper = ECIF(2016)
                    ecif = ecif_helper.get_ecif(rec_path, optimized_mol, float(6.0))
                    ld = ecif_helper.get_ligand_features_by_mol(optimized_mol)

                    data = ecif + list(ld)
                    cols = ecif_helper.get_possible_ecif() + ecif_helper.get_ligand_descriptors()
                    data_f = pd.DataFrame([data], columns=cols)
                    pred = m.predict(data_f)[0]
                    ecif_end = time.perf_counter()
                    ecif_time = ecif_time + (ecif_end - ecif_start)
                    db_helper.insert_prediction(j, optimized_mol, pred, 'results2k')
        print('following ids still cannot be processed:')
        print(err_ids)
    db_helper.__del__()
    return equibind_time, ecif_time, err_ids


def modulator(n_workers=4):
    """
    multiprocess modulator
    :param n_workers:
    :return:
    """

    only_process_err = False

    modulator_db_helper = DBUtil()
    modulator_table_name = "distinct_smiles2k"
    if not only_process_err:
        modulator_db_helper.create_err_table(modulator_table_name)
        lig_ids = modulator_db_helper.fetch_ids(modulator_table_name)
        size = math.ceil(len(lig_ids) / n_workers)
        worker_tasks = [lig_ids[i:i + size] for i in range(0, len(lig_ids), size)]
        pool = multiprocessing.Pool(n_workers)
        pool.starmap(inference_and_score_mp, zip(worker_tasks))

    # process err ids
    err_ids = modulator_db_helper.fetch_ids("err_{}".format(modulator_table_name))
    inference_and_score("err_{}".format(modulator_table_name), err_ids)


def test():
    modulator(n_workers=4)


if __name__ == '__main__':
    start = time.perf_counter()
    test()
    end = time.perf_counter()
    print('\n')
    print('run time: {} seconds'.format(round(end-start)))
