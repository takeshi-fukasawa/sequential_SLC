# Replication code for "Sequential algorithm for structural estimations with equilibrium constraints" (Fukasawa, 2026)

This repository contains replication code for the paper titled "Sequential algorithm for structural estimations with equilibrium constraints" by Takeshi Fukasawa, proposing the Sequential Linearly Constrained (SLC) algorithm for structural estimations in economics.

"mc-dynamic-game-unobs-hetero" folder contains the code for the Monte Carlo experiments on dynamic discrete games with time-varying unobserved heterogeneity. "mc_main_sequential.m" is the main code for running the numerical experiments. "mc-psd2008" folder contains the code for the Monte Carlo experiments on the convergence of the SLC/EPL algorithm and Newton's method. "check_Newton_EPL_conv.m" is the main code for running the experiments.
Note that the two replication code build on the replication code provided by Dearing and Blevins (2025). "mc-dynamic_BLP" folder contains the code for the Monte Carlo experiments on dynamic BLP (dynamic demand) estimations. "mc_sequential_main_dynamic_BLP.m" is the main code for running the experiment. The folder "sequential_algorithm" contains code for running the SLC algorithm and the Nested Fixed Point (NFXP) algorithm.
The folder "auto_diff_test" contains a document and code on the computational performance of automatic differentiation, which is briefly mentioned in the paper.

Note that the code also relies on the following repositories:
* takeshi-fukasawa/spectral (Acceleration methods of fixed-point iterations)
* takeshi-fukasawa/BLP_algorithm (Dynamic BLP)

If you run the code, please also download these repositories, and modify the path appropriately.


## References
* Dearing, A., & Blevins, J. R. (2025). Efficient and convergent sequential pseudo-likelihood estimation of dynamic discrete games. Review of Economic Studies, 92(2), 981-1021.
* Fukasawa, T. (2026). Sequential algorithm for structural estimations with equilibrium constraints. arXiv:2606.04356.
