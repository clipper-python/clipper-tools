

def print_fasta_sequence ( molecule = None ) :
    if molecule is None :
        print "Error: need to specify a molecule!"
    else :
        model = molecule.model()
        for polymer in model :
            sequence += "> chain\n"
            for monomer in polymer :
                if monomer.type().trim() == "ALA" :
                    sequence += "A"
                if monomer.type().trim() == "CYS" :
                    sequence += "C"
                if monomer.type().trim() == "ARG" :
                    sequence += "R"
                if monomer.type().trim() == "LYS" :
                    sequence += "K"
                if monomer.type().trim() == "ILE" :
                    sequence += "I"
                if monomer.type().trim() == "LEU" :
                    sequence += "L"
                if monomer.type().trim() == "ASP" :
                    sequence += "D"
                if monomer.type().trim() == "ASN" :
                    sequence += "N"
                if monomer.type().trim() == "MET" :
                    sequence += "M"
                if monomer.type().trim() == "GLY" :
                    sequence += "G"
                if monomer.type().trim() == "GLN" :
                    sequence += "Q"
                if monomer.type().trim() == "GLU" :
                    sequence += "E"
                if monomer.type().trim() == "VAL" :
                    sequence += "V"
                if monomer.type().trim() == "PHE" :
                    sequence += "F"
                if monomer.type().trim() == "TRP" :
                    sequence += "W"
                if monomer.type().trim() == "HIS" :
                    sequence += "H"
                if monomer.type().trim() == "SER" :
                    sequence += "S"
                if monomer.type().trim() == "THR" :
                    sequence += "T"
                if monomer.type().trim() == "PRO" :
                    sequence += "P"
                if monomer.type().trim() == "TYR" :
                    sequence += "Y"
                # Lots to do here!
        return sequence
