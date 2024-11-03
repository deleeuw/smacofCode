# smacofCode

Code for smacof MDS 

Ultimately this should be structured as modules. If you
want a circular configuration and your data are triads,
you combine the triad *transformation module* (currently in smacofTR) 
with the circular *Guttman module* (currently in smacofEL).

At the moment I have given up on the goal of optimizing
memory use and speed -- first I create a codebase in
R, then I'll figure out the optimization (using C).

* smacofAV - Robust MDS
* smacofBO - Bounded interval with or without additive constant
* smacofBS - Non-linear with B-splines
* smacofCR - Non-metric, Conditional rank orders
* smacofEL - Constrained, spheres and ellipsoids
* smacofHC - Rectangular, non-metric, Homogeneity Analysis, centroid restrictions
* smacofHO - Rectangular, non-metric, Homogeneity Analysis, optional rank restrictions
* smacofHS - Rectangular, non-metric, Homogeneity Analysis, spherical restrictions
* smacofID - Individual Differences
* smacofLN - Constrained, Linear Constraints
* smacofMA - Metric with additive constant
* smacofME - Metric
* smacofPC - Non-metric, Paired Comparisons
* smacofPO - Metric, Power Functions
* smacofQQ - Constrained, simplex, circumplex, uniqueness
* smacofRO - Non-metric, Rank orders
* smacofTR - Non-metric, Cartwheels
* smacofUF - Rectangular, metric, Unfolding 

