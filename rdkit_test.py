from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem import rdMolDescriptors
from chemspipy import ChemSpider


# Get molecule SMILES string from user input
smiles = input("Enter the SMILES string of the molecule: ")

# Create a molecule object from a SMILES string
mol = Chem.MolFromSmiles(smiles)

# Get Molecular Formula, Molecular Weight, and Molecular Name (InChI Key)
formula = rdMolDescriptors.CalcMolFormula(mol)
mw = Descriptors.MolWt(mol)
inchi_key = MolToInchiKey(mol)

print(f"Molecular Formula: {formula}")
print(f"Molecular Weight: {mw}")
print(f"Molecular Name (InChI Key): {inchi_key}")

# Generate 3D coordinates for the molecule
mol_3d = Chem.AddHs(mol)  # Add hydrogens to complete the valence
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # Generate 3D coordinates

# Drawing molecule
img = Draw.MolToImage(mol)  # Create an image of the molecule
img.save('example_input_output_files/molecule.png')

# Calculating Morgan Fingerprints
fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)  # Compute Morgan fingerprint with radius 2
