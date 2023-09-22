# ChebyshevHOPGD

Simulations from the manuscript Sparse grid based Chebyshev HOPGD for parameterized linear systems by Siobhán Correnty, Melina A. Freitag, and Kirk M. Soodhalter

Note, one must generate snapshots using Preconditioned Chebyshev BiCG (see https://github.com/siobhanie/ChebyshevBiCG) and then perform the tensor decomposition with HOPGD. The software for HOPGD can be obtained here: https://yelu-git.github.io/hopgd/

With this in mind, below are the scripts corresponding to the plots in the manuscript. Note, the low-rank approximations have been precomputed and saved as .mat files.

Section 4.1 2dhelm_sparse_grid.py generates the finite element matrices and plots solutions use_ChebBiCG1.m and use_ChebBiCG2.m generate Figure 4.2 See scripts make_figures_1.m and make_figures_2.m to generate Figures 4.3 and 4.4

Section 4.2 2d_helm_nonsparse.py creates the finite element matrices and plots solutions interp_nonsparse.m generates Figure 4.8b

Section 4.3 parameter_est_exper.py generates the finite element matrices parameter_est.m, parameter_est2.m, and parameter_est3.m perform the parameter estimation problem (Figure 4.12)

Section 5 interpolate_the_approx.m makes Figure 5.1 parameter_est.m, parameter_est2.m, and parameter_est3.m perform the parameter estimation (Figure 5.2)
