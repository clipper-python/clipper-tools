# Clipper-Tools
A pure Python module containing tools for crystallography and single-particle cryo-EM.
Performs highly-efficient validation metrics calculations for structure validation.

## em 
Contains functions for cutting density out of cryo-EM maps and calculating structure factors

## io
Contains functions for reading and writing maps, map coefficients, model coordinates, and structure factors

## metrics
Provides an efficient way to obtain comprehensive per-residue metrics from a model file
* **percentiles**
  * Calculate percentile rankings for per-residue metrics, compared to modern structures of similar resolution in the PDB
* **reflections**
  * Contains the *ReflectionsHandler* class, which provides a fast and simple way to get electron density at given coordinates from a reflections file
* **rotamer**
  * Employs data from the Richardson lab to calculate fit scores and discrete classifications for rotamer conformations

## utils
A number of functions that are utilised extensively within the library, and are also useful on their own, including:
* Needleman-Wunsch algorithm for pairwise sequence alignment
* Mainchain and sidechain torsion angle calculators
* B-factor analysis
* Amino acid atomic geometry validation
* Ramachandran validation
* Rotamer validation
* General calculations and matrix operations

## xray
Contains functions that make it easy to do molecular replacement with cryo-EM maps
