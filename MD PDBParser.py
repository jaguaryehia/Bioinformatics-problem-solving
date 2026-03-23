# from Bio.PDB import PDBParser
#
# def read_pdb_file(pdb_file):
#     parser = PDBParser()
#     structure = parser.get_structure('my_structure', pdb_file)
#
#     # Extract relevant information from the structure
#     coordinates = []
#     velocities = []
#     atom_types = []
#     connectivity = []
#
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 for atom in residue:
#                     # Extract coordinates
#                     coordinates.append(atom.get_coord())
#
#                     # Extract velocities (if available in the PDB file)
#                     # Note: Velocity information is not commonly provided in PDB files
#                     velocities.append(atom.get_vector().get_array())
#
#                     # Extract atom type
#                     atom_types.append(atom.get_name())
#
#                     # Extract connectivity information (e.g., bonds)
#                     # Note: PDB files may not contain explicit connectivity information
#                     #       and you may need additional methods to infer connectivity.
#                     connectivity.append(get_connectivity(atom))
#
#     return coordinates, velocities, atom_types, connectivity
#
# def get_connectivity(atom):
#     # Implement your connectivity inference method here
#     # This can be based on atomic distances, residue information, etc.
#     # Example: Assume all atoms in the same residue are connected
#     residue = atom.get_parent()
#     residue_atoms = [a.get_name() for a in residue]
#     return residue_atoms
#
# # Usage example
# pdb_file = 'pddbb/complexCID.pdb'
# coordinates, velocities, atom_types, connectivity = read_pdb_file(pdb_file)
#
# # Print the extracted information
# print('Coordinates:', coordinates)
# print('Velocities:', velocities)
# print('Atom Types:', atom_types)
# print('Connectivity:', connectivity)
# import biopandas.pdb as bpd
#
# def read_pdb_file(pdb_file):
#     pdb_data = bpd.PandasPdb().read_pdb(pdb_file)
#
#     # Extract atom information
#     atom_info = pdb_data.df['ATOM'][['atom_name', 'residue_name', 'residue_number', 'chain_id', 'x_coord', 'y_coord', 'z_coord']]
#     atom_info = atom_info.to_dict('records')
#
#     return atom_info
#
# # Example usage
# pdb_file = 'pddbb/complexCID.pdb'
# atoms = read_pdb_file(pdb_file)
#
# # Access the extracted information
# for atom in atoms:
#     print(atom['atom_name'], atom['residue_name'], atom['residue_number'], atom['chain_id'], atom['x_coord'], atom['y_coord'], atom['z_coord'])

