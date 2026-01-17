from rdkit import Chem
from rdkit.Chem import AllChem
import base64

# Define the adsorbate molecule: CH3SH (methanethiol)
smiles_CH3SH = 'CS'  # CH3SH (methanethiol)
mol_adsorbate = Chem.MolFromSmiles(smiles_CH3SH)
mol_adsorbate = Chem.AddHs(mol_adsorbate)
AllChem.EmbedMolecule(mol_adsorbate, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_adsorbate)

# Define the adsorbent surface: H-O-Zn-O ion
smiles_ion = 'O[Zn+2][O-]'
mol_adsorbent = Chem.MolFromSmiles(smiles_ion)
mol_adsorbent = Chem.AddHs(mol_adsorbent)
AllChem.EmbedMolecule(mol_adsorbent, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_adsorbent)

# Get the conformer coordinates and shift adsorbate away from adsorbent
conf_ads = mol_adsorbate.GetConformer()
for i in range(mol_adsorbate.GetNumAtoms()):
    pos = conf_ads.GetAtomPosition(i)
    conf_ads.SetAtomPosition(i, (pos.x + 4, pos.y + 2, pos.z))

# Get molecule blocks
mol_block_ads = Chem.MolToMolBlock(mol_adsorbate)
mol_block_ads_surface = Chem.MolToMolBlock(mol_adsorbent)

# Create HTML content manually to fix the zoomTo ordering issue
html_content = """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>H-O-Zn-O ion with CH3SH Adsorption</title>
    <script src="https://cdn.jsdelivr.net/npm/3dmol@2.5.3/build/3Dmol-min.js"></script>
    <style>
        body {
            margin: 0;
            padding: 10px;
            font-family: Arial, sans-serif;
        }
        #viewer {
            width: 100%;
            height: 600px;
            border: 1px solid #ccc;
            position: relative;
        }
        h1 {
            color: #333;
            margin-top: 0;
        }
    </style>
</head>
<body>
    <h1>H-O-Zn-O Ion with CH3SH Molecule Adsorption</h1>
    <div id="viewer"></div>
    
    <script>
        // Wait for 3Dmol to load
        window.addEventListener('load', function() {
            // Create viewer
            let viewer = $3Dmol.createViewer(document.getElementById('viewer'), {
                backgroundColor: 'white'
            });
            
            // Add first molecule (CH3SH adsorbate)
            let data1 = `""" + mol_block_ads.replace('\n', '\\n').replace('"', '\\"') + """`;
            viewer.addModel(data1, 'mol');
            
            // Add second molecule (H-O-Zn-O ion adsorbent)
            let data2 = `""" + mol_block_ads_surface.replace('\n', '\\n').replace('"', '\\"') + """`;
            viewer.addModel(data2, 'mol');
            
            // Set style
            viewer.setStyle({}, {sphere: {radius: 0.3}, stick: {radius: 0.15}});
            
            // Zoom to show all
            viewer.zoomTo();
            
            // Render
            viewer.render();
        });
    </script>
</body>
</html>
"""

file_name = 'H_O_Zn_O_ion_CH3SH_adsorption.html'
with open(file_name, 'w') as f:
    f.write(html_content)

print(f"Adsorption structure saved to {file_name}")
print(f"Structure shows CH3SH molecule and H-O-Zn-O ion in adsorption configuration")
