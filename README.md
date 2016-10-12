# swot (Super W Of Theta)

Author: Jean Coupon.

Contributors: Alexie Leauthaud, Martin Kilbinger, Elinor Medezinski.

## Description

**[NEW]: fits support added for input files (faster reading, names for input columns, and [CFITSIO filters](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/filters.html)).**

`SWOT` is a software written in C to compute astrophysical two-point correlation functions and weighted histograms, along with fast sky-based resampling error and covariance matrix estimations.

Available estimators include:
- angular two-point correlation function (auto and cross), "w(theta)",
- projected two-point correlation function (auto and cross), "w(R)",
- galaxy-galaxy lensing, "Delta Sigma (R)",
- 3D two-point correlation function, "xi(r)",
- 3D two-point correlation function projected along the line of sight, "xi(rp, pi), wp(rp)",
- and weighted histogram (e.g. "VMAX" estimator of mass functions).

Each estimator is wrapped around a set of "divide and conquer" algorithms, mainly (but not limited to) data storage in binary trees, approximations at large angular scale, and parellelization.

Contributions:
- the core algorithm to compute the number of pairs from a kd-tree is based on Martin Kilbinger's Athena code: http://www.cosmostat.org/software/athena/.
- the galaxy-galaxy lensing core algorithm is based on Alexie Leauthaud's code (Leauthaud et al. 2010, ApJ, 709, 97).
- Elinor Medezinski helped improving and correcting a number of bugs in the cross-correlation module.

If you use this software, please cite [Coupon et al. (2012, A&A, 542A, 5)](http://adsabs.harvard.edu/abs/2012A%26A...542A...5C). A static link is also available at http://jeancoupon.com/swot.

## Installation

### Dependencies

`SWOT` requires OPEN MPI (for paralellisation), GSL (for integration and random numbers), and CFITSIO libraires.

To install OPEN MPI, please visit http://www.open-mpi.org/. Note that OPEN MPI is a wrapper to the default C-code compiler on your machine. To install it with another C compiler (for example intel `icc`), install OPEN MPI with the following command:
```shell
$ ./configure --prefix=YOUR/PATH CC=icc CXX=icpc
```

To download and install GSL, please visit http://www.gnu.org/software/gsl/.


To download and install CFITSIO, please visit http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html.

### SWOT

1. Download the latest version here https://github.com/jcoupon/swot/releases/latest
2. Untar the archive and change directory to `swot-X.Y.Z`
3. Run `Make`

If OPEN MPI or GSL are installed in a different directory than `/usr/local`, edit the `Makefile` and set the `MPI` or `GSL` variable accordingly, or run:
```shell
$ make GSL=path/to/gsl MPI=path/to/mpi
```

The binary file is installed in swot-X.Y.Z/bin. Simply update your `PATH` variable or move `swot` to your working directory.

### VENICE

`venice` is a companion software to draw random points out of DS9 or fits masks. See https://github.com/jcoupon/venice.

## Usage


Run the software:
```
$ mpirun -np [Ncpu] swot -c configFile -corr ESTIMATOR [options]:
```

Display the default configuration file:
```
$ swot -d: display a default configuration file
```

*Important*: if using "RADEC" coordinates, the angle in the input catalogues must be in decimal degrees.

## Options


###  Estimator (`-corr` and `-est`)

- `auto`, `cross`: two-point correlation function. The choice for the estimator (`-est`) can be made among:
`ls`: Landy & Szalay (1993) estimator (default), `nat`: natural estimator: 1+DD/RR, `ham`: Hamilton (1993), or `peebles`: Peebles (1993).
- `gglens`
- `auto_wp`, `cross_wp`
- `auto_3D`, `cross_3D`
- `number`


### Open angle (`-OA`)

From Athena's documentation:
>	"If two nodes see each other under angles which are smaller than the open-angle
	threshold (OATH in config file), the tree is not further descended and those
	nodes are correlated directly. The mean, weighted galaxy ellipticities of the
	both nodes are multiplied, the angular separation is the binned barycenter
	distance.  This is of course very time-saving, which is the great advantage of a
	tree code over a brute-force approach. It introduces however errors in the
	correlation function, in particular on large scales, because the galaxy
	ellipticities are smeared out. One should play around with the value of OATH to
	see its influence on the results. As a guideline, OATH = 0.05 (about 3 degrees)
	can be used."


### Weighted quantities:

`-weighted`: compute weighted quantities ([X]: column for the weights if available). WARNING: if set, all input files must contain a weight column (including the random samples, put 1.0's if no weight).


### Input format (ascii files)


The input format is set with the `cols[1,2]` and `rancols[1,2]` options (the first column is "1"):

* default input format for `number` `RA DEC m [weight]`: change column ids with: `-cols1 1,2,3[,4]`, with `m` the quantity to compute the histogram for. Positions on the sky are required for the computation of sky-based resampling errors.

* default input format for `auto` and `cross` `[-proj como/phys] [-weighted]`: `RA DEC [z] [weight]`. This means -> `RA DEC z` (if `-proj phys`), `RA DEC w` (if `-weighted`), or `RA DEC z w` (if `-proj phys -weighted`). Change column ids with: `-cols1 1,2[,3][,4]`

* default input format for `auto_3D` and `cross_3D` `[-weighted]`: `X Y Z [weight]`. Change column ids with: `-cols1 1,2,3[,4]`. Trick: for 2D cartesian correlations, simply fix the 3rd column to a constant number.

* default input format for `auto_wp` and `cross_wp`: `RA DEC z [weight]`. Change column ids with: `-cols1 1,2,3[,4]`

* default input format for `gglens`. Lenses (`-cols1`) `RA DEC z sigz`. Change column ids with: `-cols1 1,2,3,4`. Sources (`-cols2`) `RA DEC z sigz e1 e2 weight`. Change column ids with: `-cols2 1,2,3,4,5,6,7`.

* default input format for `gglens` to compute calibration factor (`-calib yes`). Lenses (`-cols1`) `RA DEC z sigz`. Change column ids with: `-cols1 1,2,3,4`. Sources (`-cols2`) `RA DEC z sigz calib e2 weight`. Change column ids with: `-cols2 1,2,3,4,5,6,7`, where calib = 1+m or c.

### Input format (fits files)

Same as above, except that column ids may be replaced by their actual names.

`-fits`: [auto,yes,no]. If "auto", will interpret files with ".fits" or ".fit" as binary fits tables.


### bin configuration

`-nbins`: controls the number of bins.

`-nbins_pi`: controls the number of bins in the "pi" direction for wp(rp).

### tree and weights

`printTree`: outputs the tree structure in "out".data[1,2]\_tree.[ascii,fits]. Output format in ascii or fits, depending on input file format. The recorded columns are "dim1 dim2 ... dimNDIM weight1 weight2 ... weigntNsamples rank".

`printTreeAndExit`: same as above but exists after.



## Memory usage

To reduce the memory usage one can reduce the number of resamplings samples (-nsamples 8) for example, but the covariance matrix is likely to be inaccurate.

Maximum efficiency is reached when the number of cpus is a power of two.

## How to contribute to this project

(Workflow taken from https://github.com/github/gitignore)

1. [Fork this project][fork] to your account.
2. [Create a branch][branch] for the change you intend to make.
3. Make your changes to your fork.
4. [Send a pull request][pr] from your fork’s branch to our `master` branch.

Using the web-based interface to make changes is fine too, and will help you
by automatically forking the project and prompting to send a pull request too.

[fork]: http://help.github.com/forking/
[branch]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[pr]: http://help.github.com/pull-requests/
