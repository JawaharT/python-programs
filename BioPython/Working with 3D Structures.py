#.pdb or .cif files for 3D structures
from Bio.PDB import PDBParser,MMCIFParser

#Create parser
parser = PDBParser()
#Get structure from .pdb file, 1st argument is structure ID
structure = parser.get_structure("6LU7","6lu7.pdb")
print(structure)

#Models in structure
print(len(structure))

#Since there is only 1 structure
model = structure[0]

#Structure: Inside models are chains, inside chains are residue, inside residue are atoms

#Check for chains
for chain in model:
    print(f'chain {chain}, chain_ID {chain.id}')

#Check for residue
for chain in model:
    print(f'chain {chain}, chain_ID {chain.id}')
    for residue in chain:
        print(residue)

#Check for atoms
for chain in model:
    print(f'chain {chain}, chain_ID {chain.id}')
    for residue in chain:
        #print(residue)
        for atom in residue:
            print(atom)
