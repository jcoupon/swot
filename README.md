# swot (Super W Of Theta)

Authors: Jean Coupon, Alexie Leauthaud, Martin Kilbinger.

## Description

`SWOT` is a software written in C to compute astrophysical two-point correlation functions and weighted histograms, along with fast sky-based resampling error and covariance matrix estimations.

Available estimators include:
- angular two-point correlation function (auto and cross), "w(theta)",
- projected two-point correlation function (auto and cross), "w(R)",
- galaxy-galaxy lensing, "Delta Sigma (R)",
- 3D two-point correlation function, "xi(r)",
- 3D two-point correlation function projected along the line of sight, "xi(rp, pi), wp(rp)",
- and weighted histogram (e.g. "VMAX" estimator of mass functions).

Each estimator is wrapped around a set of "divide and conquer" algorithms, mainly (but not limited to) data storage in binary trees, approximations at large angular scale, and parellelization.

Contributors:
- the core algorithm to compute the number of pairs from a kd-tree is based on Martin Kilbinger's Ahtena code: http://www.cosmostat.org/software/athena/
- the galaxy-galaxy lensing core algorithm is based on Alexie Leauthaud's code (Leauthaud et al. 2010, ApJ, 709, 97).
- the rest was written by Jean Coupon.

If you use this software, please cite Coupon et al. (2012, A&A, 542A, 5).

## Installation

### Main software

1. Download the latest version here https://github.com/jcoupon/swot/releases/latest
2. Untar the archive and change directory to `swot-X.Y.Z`
3. Run `Make`

The binary file is installed in swot-X.Y.Z/bin. Simply update your `PATH` variable or move `swot` to your desired directory. 

### Dependencies

`SWOT` requires OPEN MPI (for paralellisation) and GSL (for integration and random numbers) libraires installed on your machine. If OPEN MPI or GSL is installed in a different directory than `usr/local`, edit the `Makefile` and set the `MPI` or `GSL` variable accordingly.

To install OPEN MPI, please visit http://www.open-mpi.org/. Note that OPEN MPI is a wrapper to the default C-code compiler on your machine. To install it with another C compiler (for example intel `icc`), re-install OPEN MPI with the following command:
``` 
./configure --prefix=YOUR/PATH CC=icc CXX=icpc
```

To download and install GSL, please visit http://www.gnu.org/software/gsl/

## Usage


Run the software:
```
	mpirun -np [Ncpu] swot -c configFile -corr ESTIMATOR [options]: 
```

Display the default configuration file:
```
    swot -d: display a default configuration file
```

*Important*: if using "RADEC" coordinates, the angle
in the input catalogues must be in decimal degrees.

## Options


###  Estimator
The choice for the estimator can be made among:
ls: Landy & Szalay (1993) estimator (default)
nat: Natural estimator: 1+DD/RR
ham: Hamilton (1993) estimator
peebles: Peebles (1993) estimator

* Open-angle (OA). From Athena's documentation:
"If two nodes see each other under angles which are smaller than the open-angle
threshold (OATH in config file), the tree is not further descended and those
nodes are correlated directly. The mean, weighted galaxy ellipticities of the
both nodes are multiplied, the angular separation is the binned barycenter
distance.  This is of course very time-saving, which is the great advantage of a
tree code over a brute-force approach. It introduces however errors in the
correlation function, in particular on large scales, because the galaxy
ellipticities are smeared out. One should play around with the value of OATH to
see its influence on the results. As a guideline, OATH = 0.05 (about 3 degrees)
can be used."

* default input format for "number" [-weighted]
RA DEC X [weight]
change column ids with: -cols1 1,2,3[,4]

* default input format for "auto" and "cross" [-proj como/phys] [-weighted]
RA DEC [z] [weight]
This means -> RA DEC z (if -proj phys), RA DEC w (if -weighted), 
or  RA DEC z w (if -proj phys -weighted) 
change column ids with: -cols1 1,2[,3][,4]

* default input format for "auto_3D" and "cross_3D" [-weighted]
X Y Z [weight]
change column ids with: -cols1 1,2,3[,4]
Trick: for 2D cartesian correlations, simply fix the 3rd column
to a constant number

* default input format for "auto_wp" and "cross_wp"
RA DEC z [weight]
change column ids with: -cols1 1,2,3[,4]

* default input format for "gglens"
lenses (-cols1) RA DEC z sigz
change column ids with: -cols1 1,2,3,4
sources (-cols2) RA DEC z sigz e1 e2 weight
change column ids with: -cols2 1,2,3,4,5,6,7

* default input format for "gglens" to compute calibration factor (-calib yes)
lenses (-cols1) RA DEC z sigz
change column ids with: -cols1 1,2,3,4
sources (-cols2) RA DEC z sigz calib e2 weight
change column ids with: -cols2 1,2,3,4,5,6,7
Where calib = 1+m or c

-------------------------------------------------------------------------
Compilation
-------------------------------------------------------------------------

To compile it, you need to have mpicc installed on your 
machine. Please visit http://www.open-mpi.org/ for more 
information. You also need the gsl library
(http://www.gnu.org/software/gsl/). Then run:
> make
If you want to use a different compiler than gcc (for example icc),
you need first to build mpicc with the following command:
(from intel website)
> ./configure --prefix=/usr/local CC=icc CXX=icpc F77=ifort FC=ifort
> make all install

-------------------------------------------------------------------------
Memory usage
-------------------------------------------------------------------------

To reduce the memory usage one can reduce the number of resamplings samples 
(-nsamples 8) for example, but the covariance matrix is likely to be 
inaccurate

Maximum efficiency is reached when the number of cpus is a power of two.


## How to contribute to this project

(Workflow taken from https://github.com/github/gitignore/edit/master/README.md)

1. [Fork this project][fork] to your account.
2. [Create a branch][branch] for the change you intend to make.
3. Make your changes to your fork.
4. [Send a pull request][pr] from your fork’s branch to our `master` branch.

Using the web-based interface to make changes is fine too, and will help you
by automatically forking the project and prompting to send a pull request too.

[fork]: http://help.github.com/forking/
[branch]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[pr]: http://help.github.com/pull-requests/


