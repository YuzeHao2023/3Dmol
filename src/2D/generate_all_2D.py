from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# Define all molecules/ions to generate
molecules_list = [
    ("H-O-Zn-O_ion", "[O-][Zn]O"),
    ("h2s_molecule", "S"),
    ("H_O_Zn_O_ion", "O[Zn+2][O-]"),
    ("CH3SH_molecule", "CS"),
    ("CH3_O_Zn_O_ion", "CO[Zn+2][O-]"),
    ("(CH3)2S_molecule", "CSC")
]

# Generate 2D structures for each molecule/ion
for name, smiles in molecules_list:
    print(f"Generating 2D structure for {name}...")
    
    try:
        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            print(f"  Error: Could not parse SMILES for {name}")
            continue
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Generate PNG image
        img = Draw.MolToImage(mol, size=(600, 600), includeAtomNumbers=False)
        png_filename = f"{name}_2D.png"
        img.save(png_filename)
        print(f"  ✓ PNG saved: {png_filename}")
        
        # Generate SVG image
        svg_string = Draw.MolToSVG(mol, width=600, height=600)
        svg_filename = f"{name}_2D.svg"
        with open(svg_filename, 'w') as f:
            f.write(svg_string)
        print(f"  ✓ SVG saved: {svg_filename}")
        
    except Exception as e:
        print(f"  Error generating {name}: {str(e)}")

print("\nAll 2D structures generated successfully!")
