# BME 143 - Final Model Implementation
## Samee Mushtak
An implementation in MATLAB of the mathematical model developed by Martin Howard and Pieter Rein ten Wolde in [*Finding the Center Reliably: Robust Patterns of Developmental Gene Expression* (2005)](https://doi.org/10.1103/PhysRevLett.95.208103) to examine an aspect of morphogenesis in *Drosophila* embryos.
### Function Files
* `pdefun.m` - Function defining the system of partial differential equations used in the paper's mathematical model.
* `pdefunExtended.m` - Function generalizing the system of partial differential equations defined in `pdefun.m` by allowing Bcd, corepressor, and Bcd-corepressor complex to have distinct diffusion coefficients. 
* `pdebc.m` - Function defining the boundary conditions on the embryo (zero-flux).
* `pdeic.m` - Function defining the initial densities of all compounds (zero for all).
* `normalize_density.m` - Given a matrix of protein densities, outputs a matrix with all densities normalized with respect to the maximum.
* `x_boundary.m` - Given a matrix of protein densities, calculates the position of the Hb boundary.
### Executable Files
* `pdeScript.m` - Script used to replicate figures found in paper. Solves the system of partial differential equations defined by `pdefunExtended.m`, `pdebc.m`, and `pdeic.m` 100 times to simulate 100 embryos.
* `modelExtension.m` - Script used to implement the model extension (variation of diffusion coefficients). Solves the system of partial differential equations defined by `pdefun.m`, `pdebc.m`, and `pdeic.m` 100 times to simulate 100 embryos.
* `boundaryStats.m` - Can be used to calculate the mean and standard deviation in the position of the Hb boundary for a simulation of 100 embryos using the protein density matrix generated in `pdeScript.m`
* `plots.m` - Generates the figures shown in the write-up.
