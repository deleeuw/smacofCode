# smacofCode

This repository has code in R and C for new versions of smacof. The version in the smacofRR directory is completely in R, the version in smacofRC has an R frontend and backend with the computational engine in C. The versions in the directories version_x.0 have older versions completely in C that compile to command-line executables. 

The smacof's in this repository are for a single symmetric matrix of dissimilarities (for now). They
have the following properties.

01. Transformations can be ratio, interval, polynomial, B-splinical.
02. Transformations can be either ordinal or nominal.
03. Stress and stress minimization can use weights, or not use weights.
04. Initial configurations can be read from a file, or computed by classical scaling, or random.
05. Inner knots for the B-splines can be read from a file, or equally spaced, or equally
    spaced on the percentile scale. Inner knots are simple.
06. For polynomial transformations there are no inner knots and we use the
    Bernstein basis. This has interval as a special case.
07. One outer iteration consists of a (bounded) number of inner iterations updating
    the configuration and a (bounded) number of inner iterations updating the
    transformation.
08. All matrices are stored as vectors in row-major mode. There is no double indexing.
09. Input is from a parameter file and from data files for dissimilarities and weights
    (in the weighted case) and optionally for the initial configuration and the
    inner knots of the spline. All parameters are small integers.
10. Optimization actively imposes the constraint that the configuration is centered and has unit
    sum of squares. 
11. Inner iterations for the configuration use the Guttman transform, inner iterations
    for the transformation use cyclic coordinate descent.
12. The fitted splines can be constrained to go through the origin. Combined with degree
    one this handles ratio data.
13. Monotone splines in smacof are B-splines with monotone coefficients, a.k.a. I-splines.
14. smacofRC uses a single .C() call to the computational engine. The rest is R.
15. The backend and the frontend, written in R, are shared by smacofRR and smacofRC/
