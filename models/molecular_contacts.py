

def salt_bridges ( molecule = None ) :

    if molecule is None :
        print "Error: need to specify a molecule!"
    else :
        positives = []
        negatives = []
        model = molecule.model()
        for polymer in model :
            for monomer in polymer :
                if monomer.type().trim() == "ARG" or \
                monomer.type().trim() == "LYS" :
                    positives.append (monomer)
                elif monomer.type().trim() == "ASP" or \
                monomer.type().trim() == "GLU" :
                    negatives.append(monomer)

        for positive in positives :
            for negative in negatives :
                try:
                    pos_cg = positive.find(clipper.String("CG")).coord_orth
                    neg_cg = negative.find(clipper.String("CG")).coord_orth
                    if clipper.Coord_orth.length ( pos_cg, neg_cg ) < 6.0 :
                        print positive.id(), negative.id()
                except:
                    pass
