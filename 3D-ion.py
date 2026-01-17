from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# Define the molecule using SMILES for CH3-O-Zn-O ion.
# 'CO[Zn+2][O-]' represents CH3-O-Zn-O ion with Zn in +2 oxidation state and one O with -1 charge.
smiles_string = 'CO[Zn+2][O-]'
mol = Chem.MolFromSmiles(smiles_string)

# Add hydrogens to the molecule to complete the valency
mol = Chem.AddHs(mol)

# Generate a 3D conformation for the molecule
AllChem.EmbedMolecule(mol, AllChem.ETKDG()) # Embed molecule in 3D
AllChem.MMFFOptimizeMolecule(mol) # Optimize the 3D geometry using MMFF

# Create a py3Dmol view
view = py3Dmol.view(width=400, height=400)
view.addModel(Chem.MolToMolBlock(mol), 'mol')
view.setStyle({'sphere':{'radius':0.3}, 'stick':{'radius':0.1}}) # Adjust sphere and stick radii
view.zoomTo()

# Save the view to an HTML file
html_content = view.write_html()

file_name = 'CH3_O_Zn_O_ion.html'
with open(file_name, 'w') as f:
    f.write(html_content)

print(f"3D structure saved to {file_name}")