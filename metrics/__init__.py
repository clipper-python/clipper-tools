"""
Copyright 2019 William Rochira at the University of York
Developed at York Structural Biology Laboratory - Cowtan group

- Simplifies calculation of per-residue metrics from a Clipper-Python MiniMol type
- Introduces types: MetricsModel, MetricsChain, MetricsResidue;
- MetricsModel must be initialised with a Clipper-Python MiniMol model type
- MetricsChain must be initialised with a Clipper-Python MiniMol polymer (chain) type
- MetricsResidue must be initialised with a Clipper-Python MiniMol monomer (residue) type, but also accepts arguments that contextualise it within the chain, allowing for further metric calculations, such as mainchain torsion angles and Ramachandran Plot probability
- When initialised by a MetricsChain instance, each MetricsResidue object is automatically created with all the context arguments specified
"""

import math

try:
    import clipper
    CLIPPER_MODE = 0
except ImportError:
    try:
        from clipper_python import _clipper as clipper
        CLIPPER_MODE = 1
    except ImportError:
        raise Exception('failed to import Clipper-Python')

from .. import _defs
from .. import utils
from .. import rotamer


class MetricsModel(object):
    def __init__(self, minimol):
        self._index = -1
        self.chains = [ MetricsChain(c) for c in minimol.model() ]
        self.chain_count = len(self.chains)
        self.resolution = None

    def __iter__(self):
        return self

    # Python 3: def __next__(self):
    def next(self):
        if self._index < len(self.chains)-1:
            self._index += 1
            return self.chains[self._index]
        self._index = -1
        raise StopIteration

    def add_reflections_data(self, residue_scores, mc_scores=None, sc_scores=None):
        if residue_scores is None:
            for chain in self.chains:
                for residue in chain:
                    residue.residue_fit_score = None
        else:
            for chain_id, chain_scores in enumerate(residue_scores):
                for residue_id, score in enumerate(chain_scores):
                    self.chains[chain_id].residues[residue_id].residue_fit_score = score
        if mc_scores is None:
            for chain in self.chains:
                for residue in chain:
                    residue.mainchain_fit_score = None
        else:
            for chain_id, chain_scores in enumerate(mc_scores):
                for residue_id, score in enumerate(chain_scores):
                    self.chains[chain_id].residues[residue_id].mainchain_fit_score = score
        if sc_scores is None:
            for chain in self.chains:
                for residue in chain:
                    residue.sidechain_fit_score = None
        else:
            for chain_id, chain_scores in enumerate(sc_scores):
                for residue_id, score in enumerate(chain_scores):
                    self.chains[chain_id].residues[residue_id].sidechain_fit_score = score

    def remove_chain(self, chain_id):
        matching_chains = [ chain for chain in self.chains if chain.chain_id == chain_id ]
        if len(matching_chains) == 0:
            print('Error removing chain, no chains matching that ID were found.')
        else:
            for chain in matching_chains:
                self.chains.remove(chain)
                self.chain_count -= 1


class MetricsChain(object):
    def __init__(self, mmol_chain):
        self._index = -1
        self.residues = [ ]
        self.length = len(mmol_chain)
        self.chain_id = str(mmol_chain.id().trim())
        for i, mmres in enumerate(mmol_chain):
            previous = mmol_chain[i-1] if i>0 else None
            next = mmol_chain[i+1] if i<len(mmol_chain)-1 else None
            residue = MetricsResidue(mmres, i, previous, next)
            self.residues.append(residue)
        for i, residue in enumerate(self.residues):
            if (i > 0 and i < len(self.residues)-1) and \
               (self.residues[i-1].is_aa and residue.is_aa and self.residues[i+1].is_aa) and \
               (self.residues[i-1].sequence_number+1 == residue.sequence_number == self.residues[i+1].sequence_number-1):
                residue.is_consecutive_aa = True
            else:
                residue.is_consecutive_aa = False

    def __iter__(self):
        return self

    # Python 3: def __next__(self):
    def next(self):
        if self._index < self.length-1:
            self._index += 1
            return self.residues[self._index]
        self._index = -1
        raise StopIteration

    def remove_residue(self, residue):
        if residue in self.residues:
            self.residues.remove(residue)
            self.length -= 1
        else:
            print('Error removing residue, no matching residue was found.')


class MetricsResidue(object):
    def __init__(self, mmol_residue, index_in_chain=None, previous=None, next=None):
        #self.mmol_residue = mmol_residue
        self.initialised_with_context = index_in_chain is not None
        self.index_in_chain = index_in_chain
        self.previous = previous
        self.next = next
        self.sequence_number = mmol_residue.seqnum()
        self.code = mmol_residue.type().trim()
        self.code_type = utils.code_type(mmol_residue)
        self.backbone_atoms = utils.get_backbone_atoms(mmol_residue)
        self.backbone_atoms_are_correct = None not in self.backbone_atoms
        self.backbone_geometry_is_correct = utils.check_backbone_geometry(mmol_residue) if self.backbone_atoms_are_correct else None
        self.is_aa = utils.check_is_aa(mmol_residue)
        self.is_consecutive_aa = None
        self.phi = clipper.MMonomer.protein_ramachandran_phi(self.previous, mmol_residue) if self.previous else None
        self.psi = clipper.MMonomer.protein_ramachandran_psi(mmol_residue, self.next) if self.next else None
        if self.phi is not None and math.isnan(self.phi):
            self.phi = None
        if self.psi is not None and math.isnan(self.psi):
            self.psi = None
        self.chis = utils.calculate_chis(mmol_residue)
        self.is_sidechain_complete = _defs.SC_INCOMPLETE_STRING not in self.chis
        self.ramachandran_score = utils.calculate_ramachandran_score(mmol_residue, self.phi, self.psi)
        self.rotamer_score = utils.calculate_rotamer_score(mmol_residue, chis=self.chis) if self.is_sidechain_complete else None
        self.max_b_factor, self.avg_b_factor, self.std_b_factor = utils.analyse_b_factors(mmol_residue)
        self.residue_fit_score = None
        self.mainchain_fit_score = None
        self.sidechain_fit_score = None
