from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# Define the adsorbate molecule: CH3SH (methanethiol)
smiles_CH3SH = 'CS'
mol_adsorbate = Chem.MolFromSmiles(smiles_CH3SH)
mol_adsorbate = Chem.AddHs(mol_adsorbate)
AllChem.Compute2DCoords(mol_adsorbate)

# Define the adsorbent ion: H-O-Zn-O
smiles_ion = 'O[Zn+2][O-]'
mol_adsorbent = Chem.MolFromSmiles(smiles_ion)
mol_adsorbent = Chem.AddHs(mol_adsorbent)
AllChem.Compute2DCoords(mol_adsorbent)

# Generate 2D images for both molecules
img_adsorbate = Draw.MolToImage(mol_adsorbate, size=(400, 400), includeAtomNumbers=False)
img_adsorbent = Draw.MolToImage(mol_adsorbent, size=(400, 400), includeAtomNumbers=False)

# Generate side-by-side comparison image
from PIL import Image
img_combined = Image.new('RGB', (850, 450), color='white')
img_combined.paste(img_adsorbate, (10, 25))
img_combined.paste(img_adsorbent, (440, 25))

# Add labels to identify the molecules
from PIL import ImageDraw, ImageFont
draw = ImageDraw.Draw(img_combined)
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 16)
except:
    font = ImageFont.load_default()

draw.text((150, 10), "CH3SH", fill='black', font=font)
draw.text((550, 10), "H-O-Zn-O Ion", fill='black', font=font)

# Save the combined image
png_file_name = '2D-adsorption_combined.png'
img_combined.save(png_file_name)
print(f"2D adsorption structures saved to {png_file_name}")

# Also generate individual SVG files
svg_adsorbate = Draw.MolToSVG(mol_adsorbate, width=400, height=400)
svg_adsorbent = Draw.MolToSVG(mol_adsorbent, width=400, height=400)

with open('2D-adsorption_CH3SH.svg', 'w') as f:
    f.write(svg_adsorbate)
print(f"2D CH3SH structure SVG saved to 2D-adsorption_CH3SH.svg")

with open('2D-adsorption_ion.svg', 'w') as f:
    f.write(svg_adsorbent)
print(f"2D H-O-Zn-O ion structure SVG saved to 2D-adsorption_ion.svg")

# Generate single image with both molecules
molecules = [mol_adsorbate, mol_adsorbent]
legends = ["CH3SH\n(Adsorbate)", "H-O-Zn-O Ion\n(Adsorbent)"]
img_grid = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(400, 400), legends=legends)
img_grid.save('2D-adsorption_structures.png')
print(f"2D adsorption grid image saved to 2D-adsorption_structures.png")
