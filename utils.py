"""
Copyright 2019 William Rochira at the University of York
Developed at York Structural Biology Laboratory - Cowtan group
"""
import math

from clipper_python import _clipper as clipper

from . import _defs
from . import rotamer


THREE_LETTER_CODES = { 0 : [ 'ALA', 'GLY', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP', 'SER',
                             'THR', 'CYS', 'MET', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU' ],
                       1 : [ 'MSE', 'SEC' ],
                       2 : [ 'UNK' ] }

CHI_ATOMS = [ { ('N', 'CA', 'CB', 'CG') : ('ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'TRP', 'TYR', 'MSE'),
                ('N', 'CA', 'CB', 'CG1') : ('ILE', 'VAL'),
                ('N', 'CA', 'CB', 'SG') : ('CYS'),
                ('N', 'CA', 'CB', 'SE') : ('SEC'),
                ('N', 'CA', 'CB', 'OG') : ('SER'),
                ('N', 'CA', 'CB', 'OG1') : ('THR') },
              { ('CA', 'CB', 'CG', 'CD') : ('ARG', 'GLN', 'GLU', 'LYS', 'PRO'),
                ('CA', 'CB', 'CG', 'CD1') : ('LEU', 'PHE', 'TRP', 'TYR'),
                ('CA', 'CB', 'CG', 'OD1') : ('ASN', 'ASP'),
                ('CA', 'CB', 'CG', 'ND1') : ('HIS'),
                ('CA', 'CB', 'CG1', 'CD1') : ('ILE'),
                ('CA', 'CB', 'CG', 'SD') : ('MET'),
                ('CA', 'CB', 'CG', 'SE') : ('MSE') },
              { ('CB', 'CG', 'CD', 'OE1') : ('GLN', 'GLU'),
                ('CB', 'CG', 'CD', 'NE') : ('ARG'),
                ('CB', 'CG', 'CD', 'CE') : ('LYS'),
                ('CB', 'CG', 'SD', 'CE') : ('MET'),
                ('CB', 'CG', 'SE', 'CE') : ('MSE') },
              { ('CG', 'CD', 'NE', 'CZ') : ('ARG'),
                ('CG', 'CD', 'CE', 'NZ') : ('LYS') },
              { ('CD', 'NE', 'CZ', 'NH1') : ('ARG') } ]


# General calculations
def mean(values):
    return float(sum(values)) / len(values)

def median(values):
    values = sorted(values)
    n = len(values)
    if n < 1:
        return None
    i = n//2
    if n % 2 == 1:
        return values[i]
    return mean(values[i-1:i+1])

# Matrix operations
def dot_product(xyz1, xyz2):
    return xyz1[0] * xyz2[0] + xyz1[1] * xyz2[1] + xyz1[2] * xyz2[2]

def cross_product(xyz1, xyz2):
    return [ xyz1[1] * xyz2[2] - xyz1[2] * xyz2[1],
             xyz1[2] * xyz2[0] - xyz1[0] * xyz2[2],
             xyz1[0] * xyz2[1] - xyz1[1] * xyz2[0] ]

def magnitude(xyz):
    return (xyz[0]**2 + xyz[1]**2 + xyz[2]**2) ** 0.5

def unit(xyz):
    length = magnitude(xyz)
    return [ xyz[0] / length, xyz[1] / length, xyz[2] / length ]

def subtract(xyz1, xyz2):
    return [ xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2] ]

def distance(xyz1, xyz2):
    v = subtract(xyz1, xyz2)
    return magnitude(v)

def angle(xyz1, xyz2, xyz3):
    v1 = subtract(xyz2, xyz1)
    v2 = subtract(xyz2, xyz3)
    angle = math.acos(dot_product(v1, v2) / (magnitude(v1) * magnitude(v2)))
    return math.degrees(angle)

def torsion(xyz1, xyz2, xyz3, xyz4, range_positive=True):
    b1 = subtract(xyz2, xyz1)
    b2 = subtract(xyz3, xyz2)
    b3 = subtract(xyz4, xyz3)
    n1 = cross_product(b1, b2)
    n2 = cross_product(b2, b3)
    m1 = cross_product(n1, n2)
    y = dot_product(m1, unit(b2))
    x = dot_product(n1, n2)
    angle = math.degrees(math.atan2(y, x))
    if range_positive and angle < 0:
        angle += 360
    elif not range_positive and angle > 180:
        angle -= 360
    return angle

# (MiniMol) residue functions
def code_type(mmol_residue):
    try:
        return next(category for category, group in THREE_LETTER_CODES.items() if mmol_residue.type().trim() in group)
    except:
        return None

def get_backbone_atoms(mmol_residue):
    try:
        n = next(atom for atom in mmol_residue if atom.id().trim().replace(' ', '') == 'N' or atom.id().trim().replace(' ', '') == 'N:A')
    except:
        n = None
    try:
        ca = next(atom for atom in mmol_residue if atom.id().trim().replace(' ', '') == 'CA' or atom.id().trim().replace(' ', '') == 'CA:A')
    except:
        ca = None
    try:
        c = next(atom for atom in mmol_residue if atom.id().trim().replace(' ', '') == 'C' or atom.id().trim().replace(' ', '') == 'C:A')
    except:
        c = None
    return n, ca, c

def check_backbone_geometry(mmol_residue):
    n, ca, c = get_backbone_atoms(mmol_residue)
    return distance(n.coord, ca.coord) < 1.8 and distance(ca.coord, c.coord) < 1.8

def calculate_chis(mmol_residue):
    chis = [ ]
    for i in range(5):
        chi_atoms = [ ]
        has_chi = any(mmol_residue.type().trim() in residues for residues in list(CHI_ATOMS[i].values()))
        if not has_chi:
            chis.append(None)
            continue
        required_atom_names = next(atoms for atoms, residues in CHI_ATOMS[i].items() if mmol_residue.type().trim() in residues)
        missing_atom_names = [ ]
        for required_atom_name in required_atom_names:
            found = False
            for atom in mmol_residue:
                atom_name = atom.id().trim().replace(' ', '')
                if atom_name in (required_atom_name, required_atom_name + ':A'):
                    chi_atoms.append(atom)
                    found = True
            if not found:
                missing_atom_names.append(required_atom_name)
        if len(chi_atoms) < 4:
            chis.append(_defs.SC_INCOMPLETE_STRING)
            continue
        xyzs = [ atom.coord for atom in chi_atoms ]
        chis.append(torsion(xyzs[0], xyzs[1], xyzs[2], xyzs[3]))
    return tuple(chis)

def analyse_b_factors(mmol_residue):
    residue_b_factors = [ clipper.Util_u2b(atom.u_iso) for atom in mmol_residue ]
    b_max = max(residue_b_factors)
    b_average = mean(residue_b_factors)
    b_stdev = (sum([ (x - b_average) ** 2 for x in residue_b_factors ]) / len(residue_b_factors)) ** 0.5
    return b_max, b_average, b_stdev

def check_is_aa(mmol_residue, strict=False):
    allowed_types = (0,) if strict else (0, 1)
    if code_type(mmol_residue) in allowed_types and \
       None not in get_backbone_atoms(mmol_residue) and \
       check_backbone_geometry(mmol_residue):
        return True
    return False

def calculate_rotamer_probability(mmol_residue, chis=None):
    if chis is None:
        chis = calculate_chis(mmol_residue)
    return rotamer.get_probability(mmol_residue.type().trim(), chis)

def calculate_ramachandran_probability(mmol_residue, phi, psi):
    if phi is None or psi is None:
        return None
    rama_function = clipper.Ramachandran(clipper.Ramachandran.Gly5) if mmol_residue.type().trim() == 'GLY' else \
                    clipper.Ramachandran(clipper.Ramachandran.Pro5) if mmol_residue.type().trim() == 'PRO' else \
                    clipper.Ramachandran(clipper.Ramachandran.NonGlyPro5)
    return rama_function.probability(phi, psi)
