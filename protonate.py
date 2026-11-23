from rdkit import Chem
from rdkit.Chem import AllChem
import sys

infile, outfile = sys.argv[1], sys.argv[2]
mol = Chem.SDMolSupplier(infile, removeHs=False)[0]
print("input atoms:", mol.GetNumAtoms(), Chem.rdMolDescriptors.CalcMolFormula(mol))

# strip H, protonate the primary amine, put all H back with 3D coordinates
mol = Chem.RemoveHs(mol)
matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3;H2;!$(N-C=O)]"))
assert len(matches) == 1, f"expected exactly 1 primary amine, found {len(matches)}"
n = mol.GetAtomWithIdx(matches[0][0])
n.SetFormalCharge(1)
n.SetNumExplicitHs(0)
n.SetNoImplicit(False)
Chem.SanitizeMol(mol)
mol = Chem.AddHs(mol, addCoords=True)

print("output atoms:", mol.GetNumAtoms(), Chem.rdMolDescriptors.CalcMolFormula(mol))
print("net charge:", Chem.GetFormalCharge(mol))
Chem.MolToMolFile(mol, outfile)
