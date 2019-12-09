# Clipper-Tools
A pure Python module containing simple tools for crystallography and single-particle cryo-EM

## em 
...

## io
...

## metrics
- Simplifies calculation of per-residue metrics from a Clipper-Python MiniMol type
- Introduces types: MetricsModel, MetricsChain, MetricsResidue;
- MetricsModel must be initialised with a Clipper-Python MiniMol model type
- MetricsChain must be initialised with a Clipper-Python MiniMol polymer (chain) type
- MetricsResidue must be initialised with a Clipper-Python MiniMol monomer (residue) type, but also accepts arguments that contextualise it within the chain, allowing for further metric calculations, such as mainchain torsion angles and Ramachandran Plot probability
- When initialised by a MetricsChain instance, each MetricsResidue object is automatically created with all the context arguments specified

## rotamer
...

## utils
...

## xray
...
