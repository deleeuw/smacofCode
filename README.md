# smacofCode

This repository has code in R and C for new versions of smacof. The version in the smacofRR directory is completely in R, the version in smacofRC has an R frontend and backend with the computational engine in C. The versions in the directories version_x.0 have older versions completely in C that compile to command-line executables. 

The smacof's in this repository are for a single symmetric matrix of dissimilarities (for now). They
have the following properties.

01. Transformations can be ratio, interval, polynomial, B-splinical.
02. Transformations can be either ordinal or nominal.
03. Stress and stress minimization can use weights, or not use weights.
04. Initial configurations can be given, or torgerson, or random.
05. Inner knots for the B-splines can be given, or equally spaced, or equally
    spaced on the percentile scale.
06. For polynomial transformations there are no inner knots and we use the
    Bernstein basis.
07. One outer iteration consists of a (bounded) number of inner iterations updating
    the configuration and a (bounded) number of inner iterations updating the
    transformation.
08. All matrices are stored as vectors in row-major mode. There is no double indexing.
09. Input is from a parameter file and from data files for dissimilarties and weights
    (in the weighted case) and optionally for the initial configuration and the
    inner knots of the spline.
10. Optimization imposes the constraint that the configuration is centered and has unit
    sum of squares. 
11. Inner iterations for the configuration use the Guttman transform, inner iterations
    for the transformation use cyclic coordinate descent.
12. The fitted splines can be constrained to go through the origin.
13. Monotone splines in smacof are B-splines with monotone coefficients, a.k.a. I-splines.
14. smacofRC uses a single .C() call to the computational engine. The rest is R.
