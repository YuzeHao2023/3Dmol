from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# Define the molecule using SMILES for H-O-Zn-O- ion.
# '[O-][Zn]O' represents Zinc bonded to an oxygen with a negative charge and another neutral oxygen.
# After adding hydrogens, it will form 'HO-[Zn]-O-'.
smiles_string = '[O-][Zn]O'
mol = Chem.MolFromSmiles(smiles_string)

# Add hydrogens to the molecule to complete the valency of the neutral oxygen
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

file_name = 'H-O-Zn-O_ion.html'
with open(file_name, 'w') as f:
    f.write(html_content)

print(f"3D structure saved to {file_name}")