from clipper_python import _clipper as clipper

## convenience functions ##
## please add more stuff ##

def is_aminoacid ( residue_name = "None" ) :
    residue_names={ 'UNK', 'ALA', 'GLY', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'MET', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU' }

    if residue_name == "None" :
        return False

    if residue_name in residue_names :
        return True
    else:
        return False


def is_mainchain ( atom_name = "None" ) :

    atom_names = { 'C', 'O', 'CA', 'N' }

    if atom_name in atom_names :
        return True
    else :
        return False
