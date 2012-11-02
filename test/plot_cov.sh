#! /bin/sh

MIN=0.001
MAX=1.5
FILEIN=$1
NBINS=`wc -l $FILEIN | awk '{print $1}'`
NAME=$2

#Last option: 1: coeff, 0: covariance

idl << EOF

.r plot_cov.pro
plot_cov, 1, '$FILEIN', $NBINS, $MIN, $MAX, '$NAME',1
EOF

gv graph.ps &
