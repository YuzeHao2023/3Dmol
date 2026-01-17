from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from PIL import Image, ImageDraw, ImageFont

# Define the ion using SMILES with charge information
# 修改这个 SMILES 字符串来生成不同的离子
smiles_string = 'O[Zn+2][O-]'
mol = Chem.MolFromSmiles(smiles_string)

# Add hydrogens to the molecule to complete the valency
mol = Chem.AddHs(mol)

# Generate 2D coordinates for the ion
AllChem.Compute2DCoords(mol)

# Generate 2D image of the ion
img = Draw.MolToImage(mol, size=(600, 600), includeAtomNumbers=False)

# Save the 2D structure as PNG
png_file_name = '3D-ion_2D.png'
img.save(png_file_name)

print(f"2D ion structure saved to {png_file_name}")

# Also generate SVG format for better quality
svg_string = Draw.MolToSVG(mol, width=600, height=600)
svg_file_name = '3D-ion_2D.svg'
with open(svg_file_name, 'w') as f:
    f.write(svg_string)

print(f"2D ion structure SVG saved to {svg_file_name}")
