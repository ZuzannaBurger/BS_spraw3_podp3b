import argparse
from Bio import PDB

def parse_arguments():
    p = argparse.ArgumentParser()
    p.add_argument("-i", required=True)
    p.add_argument("-o", required=True)
    return p.parse_args()

def load_templates():
    pdb_parser = PDB.PDBParser()
    return {
        "A": pdb_parser.get_structure("A", "templates/adenine.pdb")[0]["A"][26],
        "U": pdb_parser.get_structure("U", "templates/uracil.pdb")[0]["A"][11],
        "G": pdb_parser.get_structure("G", "templates/guanine.pdb")[0]["A"][24],
        "C": pdb_parser.get_structure("C", "templates/cytosine.pdb")[0]["A"][29],
    }

def reconstruct_structure(input_path, templates):
    pdb_parser = PDB.PDBParser()
    coarse_structure = pdb_parser.get_structure("430d_cg", input_path)
    full_structure = PDB.Structure.Structure("430d_cg")
    superimposer = PDB.Superimposer()

    for model in coarse_structure:
        new_model = PDB.Model.Model(model.id)
        full_structure.add(new_model)
        for chain in model:
            new_chain = PDB.Chain.Chain(chain.id)
            new_model.add(new_chain)
            for i, residue in enumerate(chain):
                res_name = residue.get_resname()
                if res_name not in templates:
                    continue
                template_residue = templates[res_name].copy()
                template_residue.id = (" ", i + 1, " ")
                target_atoms = [residue[atom.get_name()] for atom in residue]
                template_atoms = [template_residue[atom.get_name()] for atom in residue]
                superimposer.set_atoms(target_atoms, template_atoms)
                for atom in template_residue:
                    atom.transform(superimposer.rotran[0], superimposer.rotran[1])
                new_chain.add(template_residue)

    return full_structure

def save_structure(structure, output_path):
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_path)

args = parse_arguments()
templates = load_templates()
full_structure = reconstruct_structure(args.i, templates)
save_structure(full_structure, args.o)

