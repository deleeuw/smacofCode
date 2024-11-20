# smacofCode

Code for smacof MDS 

Ultimately this should be structured as modules. If you
want a circular configuration and your data are triads,
you combine the triad *transformation module* (currently in smacofTR) 
with the circular *Guttmnan module* (currently in smacofEL).

At the moment I have given up on the goal of optimizing
memory use and speed -- first I create a codebase in
R, then I'll figure out the optimization (using C).

* smacofAC - Metric with additive constant
* smacofAV - Robust MDS
* smacofBA - Smacof Basics (Introduction)
* smacofBO - Metric with bounds
* smacofBS - Non-linear with B-splines
* smacofCR - Non-metric, Conditional Rank Orders
* smacofDT - Data format
* smacofEL - Constrained, spheres and ellipsoids
* smacofFA - Accelerated Smacof
* smacofHC - Rectangular, non-metric, Homogeneity Analysis, Centroid Constraints
* smacofHO - Rectangular, non-metric, Homogeneity Analysis, Rank Constraints
* smacofHS - Rectangular, non-metric, Homogeneity Analysis, Circular Constraints
* smacofID - Individual Differences
* smacofLN - Constrained, Linear Constraints
* smacofME - Metric, with acceleration
* smacofPC - Non-metric, Paired Comparisions
* smacofPO - Non-linear, Power Functions
* smacofQQ - Constrained, simplex, circumplex, uniqueness
* smacofRO - Non-metric, Rank orders
* smacofSE - Initial, Elegant, Strain
* smacofTR - Non-metric, Triads and Propellors
* smacofUF - Rectangular, metric, Unfolding 
* smacofUT - Various Utilities
