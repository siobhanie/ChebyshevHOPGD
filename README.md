# ChebyshevHOPGD

Simulations from the manuscript Sparse grid based Chebyshev HOPGD for parameterized linear systems by Siobhán Correnty, Melina A. Freitag, and Kirk M. Soodhalter

One can generate snapshots using Preconditioned Chebyshev BiCG (see also https://github.com/siobhanie/ChebyshevBiCG) and then perform the tensor decomposition with HOPGD. The snapshots for the simulations in this work are precomputed and saved in the corresponding directories. The tensor decompositions have also been precomputed with HOPGD and saved. The software for HOPGD was not written by the authors of this paper, but rather the original developers (see: https://yelu-git.github.io/hopgd/).  

With this in mind, below are the scripts corresponding to the plots in the manuscript. In most cases, the finite element matrices have been precomputed with FEniCS and saved. The exception is the finite element matrices in the simulations in section 4.2, which must be computed and saved by the user. 

Simulations section 4.1: 2dhelm_sparse_grid.py generates the finite element matrices and plots solutions, use_ChebBiCG1.m and use_ChebBiCG2.m generate Figure 4.2. See scripts make_figures_1.m and make_figures_2.m to generate Figures 4.3 and 4.4.

Simulations section 4.2: 2d_helm_nonsparse.py creates the finite element matrices and plots solutions interp_nonsparse.m generates Figure 4.8b.

Simulations section 4.3: parameter_est_exper.py generates the finite element matrices, and the files parameter_est.m, parameter_est2.m, and parameter_est3.m perform the parameter estimation problem (Figure 4.12)

Simulations section 5: interpolate_the_approx.m makes Figure 5.1 and the set parameter_est.m, parameter_est2.m, and parameter_est3.m perform the parameter estimation (Figure 5.2)
