# smacofCode

Code for smacof MDS 

Ultimately this should be structured as modules. If you
want a circular configuration and your data are triads,
you combine the triad *transformation module* (currently in smacofTR) 
with the circular *Guttmnan module* (currently in smacofEL).

At the moment I have given up on the goal of optimizing
memory use and speed -- first I create a codebase in
R, then I'll figure out the optimization (using C).

* smacofAC - Metric with bounds and/or additive constant
* smacofBS - Non-linear with B-splines
* smacofCR - Non-metric, Conditional rank orders
* smacofEL - Constrained, spheres and ellipsoids
* smacofHO - Rectangular, non-metric, Homogeneity Analysis
* smacofID - Individual Differences
* smacofLN - Constrained, Linear Constraints
* smacofPC - Non-metric, Paired Comparisons
* smacofQQ - Constrained, simplex, circumplex, uniqueness
* smacofRO - Non-metric, Rank orders
* smacofTR - Non-metric, Triads and Propellors
* smacofUF - Rectangular, metric, Unfolding 

