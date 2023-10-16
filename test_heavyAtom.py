from rdkit import Chem
import rdkit.Chem.Descriptors as Desc

class Test:
    def __init__(self, smiles: str):
        self.smiles = smiles

    def get_heavy_atom_count(self):
        mol = Chem.MolFromSmiles(self.smiles)
        return Desc.HeavyAtomCount(mol)

if __name__ == '__main__':
    test_instance = Test('CCO')
    print(test_instance.get_heavy_atom_count())
