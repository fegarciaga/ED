# ED
This is a code for exact diagonalization of the t-t' hubbard model in up to 3 dimensions. I developed this code to compare exact results with approximate results from the constrained path quantum Monte Carlo method.

The main file y H_exact, where one set the number of fermions, the geometry of the lattice, up to three nearest-neighbor hopping parameters, up to three next-nearest-neighbor hopping parameter and up to three angle parameters for the twist angle boundary condition. The code is general and detects the dimensionality in order to set next-nearest-neighbor in a plane or an axis.

Due to non inherent non polynomial computation time, the code performs quickly on for up to 10 sites.
