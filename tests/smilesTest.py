from rdkit import Chem
from rdkit.Chem import AllChem
smiles = 'CC(C)CC(NC(=O)C(Cc1ccccc1)NC(=O)C(N)CO)C(=O)NC(CC(C)C)C(=O)NC(CCCNC(=N)N)C(=O)NC(CC(N)=O)C(=O)O'


def test():
    mol = AllChem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=False, useBasicKnowledge=False)
    conf = mol.GetConformer()
    print(conf)


if __name__ == '__main__':
    test()
