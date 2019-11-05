import numpy as np
import coordinate_kicks
import molecular_contacts
import sequences
import ramachandran

## convenience functions for models  ##
## please add more stuff as required ##

def trim_to_c_alphas ( molecule = None ) :
    if molecule is None :
        print "Error: need to specify a molecule!"
    else :
        model = molecule.model()
        new_molecule = clipper.MiniMol()
        new_model = new_molecule.model()
        for polymer in model :
            new_polymer = clipper.MPolymer()
            for monomer in polymer :
                new_monomer = clipper.MMonomer()
                for atom in monomer :
                    if "CA" in atom.name() :
                        new_monomer.insert(atom)
                        new_monomer.set_id(monomer.id())
                        new_monomer.set_seqnum(monomer.seqnum())
                        new_monomer.set_type(monomer.type())
                new_polymer.insert(new_monomer)
            new_model.insert(new_polymer)
    return new_molecule


def b_averages ( molecule = None ) :
    if molecule is None :
        print "Error: need to specify a molecule!"
    else :
        model = molecule.model()
        list_b_wat = []
        list_b_rst = []

        for polymer in model :
            for monomer in polymer :
                for atom in monomer :
                    if monomer.type().trim() == "HOH" :
                        list_b_wat.append(clipper.Util.u2b(atom.u_iso))
                    else :
                        list_b_rst.append(clipper.Util.u2b(atom.u_iso))

        if len(list_b_wat) == 0 :
          list_b_wat.append( 0.0)
        print "<B>(std) protein and ligands"
        print str(np.mean(list_b_rst)), str(np.std(list_b_rst))
        print "<B>(std) waters"
        print str(np.mean(list_b_wat)), str(np.std(list_b_wat))
