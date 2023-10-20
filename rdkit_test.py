import os
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem import rdMolDescriptors

# Get the API key securely
api_key = ''
if not api_key:
    raise ValueError("Please set the OPENAI_API_KEY environment variable.")

# Corrected endpoint for chat models
endpoint = 'https://api.openai.com/v1/chat/completions'

# Get molecule SMILES string from user input
smiles = input("Enter the SMILES string of the molecule: ")

# Create a molecule object from a SMILES string
mol = Chem.MolFromSmiles(smiles)

# Get Molecular Formula, Molecular Weight, and Molecular Name (InChI Key)
formula = rdMolDescriptors.CalcMolFormula(mol)
mw = Descriptors.MolWt(mol)
inchi_key = MolToInchiKey(mol)

# Print out the molecular details
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

# Function to communicate with ChatGPT API
def get_chatgpt_response(prompt):
    headers = {
        'Authorization': f'Bearer {api_key}',
        'Content-Type': 'application/json',
    }
    data = {
        'model': 'gpt-3.5-turbo',
        'messages': [{'role': 'user', 'content': prompt}]
    }
    response = requests.post(endpoint, headers=headers, json=data)
    if response.status_code == 200:
        return response.json()['choices'][0]['message']['content']
    else:
        error_msg = f'Error {response.status_code}: {response.json().get("error", {}).get("message", "Unknown error")}'
        print(error_msg)
        return None

# Ask the user for an additional question about the molecule
user_question = input("Please enter your question about the molecule: ")
response_text = get_chatgpt_response(user_question)
print(response_text)
