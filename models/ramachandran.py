

def ramachandran_maps ( molecule = None ) :

    rama_gly = clipper.Ramachandran ( clipper.Ramachandran.Gly5 )
    rama_pro = clipper.Ramachandran ( clipper.Ramachandran.Pro5 )
    rama_rest= clipper.Ramachandran ( clipper.Ramachandran.NonGlyPro5 )

    rama_threshold_preferred = 0.02
    rama_threshold_allowed = 0.002

    rama_gly.set_thresholds ( rama_threshold_preferred, rama_threshold_allowed )
    rama_pro.set_thresholds ( rama_threshold_preferred, rama_threshold_allowed )
    rama_rest.set_thresholds ( rama_threshold_preferred, rama_threshold_allowed )

    prev_residue = clipper.MMonomer()
    n_residues = n_residues_chain = n_allowed = n_favoured = n_outliers = 0
    phi = psi = 0.0

    model = molecule.model()

    for chain in model :
        n_residues_chain = 0
        for residue in chain :
            n_residues_chain += 1

            if n_residues_chain == 1 :
                # first residue initialises prev_residue
                prev_residue = residue
            elif n_residues_chain == len ( chain ) :
                break
            elif is_aminoacid ( prev_residue.type().trim() ) and \
                 is_aminoacid ( residue.type().trim() ) and \
                 is_aminoacid (chain [ n_residues_chain ].type().trim()) and \
                 prev_residue.seqnum() +2 == residue.seqnum() +1 == chain [ n_residues_chain ].seqnum() :
                 # need to check that they're also consecutive
                phi = clipper.MMonomer.protein_ramachandran_phi ( prev_residue, residue )
                psi = clipper.MMonomer.protein_ramachandran_psi ( residue, chain [ n_residues_chain ] )
                n_residues += 1

                if residue.type().trim() == "GLY" :
                    if rama_gly.favored ( phi, psi ) :
                        print "Favored: ", residue.type(), residue.id(), phi, psi
                        n_favoured += 1
                    elif rama_gly.allowed ( phi, psi ) :
                        print "Allowed: ", residue.type(), residue.id(), phi, psi
                        n_allowed += 1
                    else :
                        print "Outlier: ", residue.type(), residue.id(), phi, psi
                        n_outliers += 1

                elif residue.type().trim() == "PRO" :
                    if rama_pro.favored ( phi, psi ) :
                        print "Favored: ", residue.type(), residue.id(), phi, psi
                        n_favoured += 1
                    elif rama_gly.allowed ( phi, psi ) :
                        print "Allowed: ", residue.type(), residue.id(), phi, psi
                        n_allowed += 1
                    else :
                        print "Outlier: ", residue.type(), residue.id(), phi, psi
                        n_outliers += 1
                else :
                    if rama_rest.favored ( phi, psi ) :
                        print "Favored: ", residue.type(), residue.id(), phi, psi
                        n_favoured += 1
                    elif rama_rest.allowed ( phi, psi ) :
                        print "Allowed: ", residue.type(), residue.id(), phi, psi
                        n_allowed += 1
                    else :
                        print "Outlier: ", residue.type(), residue.id(), phi, psi
                        n_outliers += 1

            elif is_aminoacid ( prev_residue.type().trim() ) and \
                 is_aminoacid ( residue.type().trim() ) and \
                 is_aminoacid ( chain [ n_residues_chain ].type().trim()) :
                # make plans for supporting insertion codes?
                print "REJECTED triplet\n"

            prev_residue = residue

    print "Number of favored: ", n_favoured
    print "Number of allowed: ", n_allowed
    print "Number of outliers: ", n_outliers


if len(sys.argv) == 3 :
    molecule = read_molecule()
    print_fasta_sequence ( molecule = molecule )
    new_molecule = trim_to_c_alphas ( molecule = molecule )
    b_averages (molecule)
    salt_bridges (molecule)
    ramachandran_maps (molecule)
    write_molecule ( molecule = new_molecule )
else :
    print "Error: the script takes two parameters: an input file, and an output file."
    print "         example: python script_2.py input.pdb output.pdb "
