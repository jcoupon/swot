#! /bin/bash
set -o nounset

# default mpirun and swot locations
MPIRUN=mpirun
SWOT=$HOME/local/source/GitHub/swot/bin/swot

# default number of cores
NP=2

if [[ "$HOST" =~ "ojingo" ]]; then
   MPIRUN=/opt/openmpi-1.8.6_clang/bin/mpirun
   NP=8
fi
if [[ "$HOST" =~ "hayabusa" ]]; then
   MPIRUN=/opt/openmpi-1.8.5/bin/mpirun
   NP=4
fi

# ...  add your machine, mpi and swot path

main() {

   mkdir -p results

   SWOT_DEFAULT=results/swot_default.ascii

   # create default parameter file
   echo "Writing the default configuration file in $SWOT_DEFAULT ..."
   $SWOT -d > $SWOT_DEFAULT

   # run the angular two-point correlation function
   echo "Computing the angular two-point correlation function..."
   runWtheta

   # runWp

   return
}

# ----------------------------------------------------------------------------------- #
# functions
# ----------------------------------------------------------------------------------- #

function runWtheta
{

   $MPIRUN -np $NP $SWOT -c $SWOT_DEFAULT -data1 data/data1.ascii -ran1 data/ran1.ascii -o results/wtheta.out

   return
}

function runWp
{

   RAN_FILE=data/ran1.ascii
   if false; then

	   N=$( awk '!(/^#/) {n+=1} END {print n}' $RAN_FILE )

	   venice -r -xmin 0 -xmax 1 -ymin 0 -ymax 1 -coord spher -z 0.5,1.0 -npart $N \
         | awk 'NR > 1 {print $3}' | paste $RAN_FILE - \
	      | awk '!(/^#/) && $5 == 0 && $12 == 0 {print $1, $2, $15, 1.0}' \
	        > ${RAN_FILE%.ascii}_wp.ascii

   fi

   $MPIRUN -np $NP $SWOT -c $SWOT_DEFAULT -corr auto_wp -weighted  \
      -cols1 1,2,3,7 -rancols1 1,2,3,4 -range 0.2,10.0 -nbins 10 -data1 data/data1_wp.ascii -ran1 ${RAN_FILE%.ascii}_wp.ascii -proj phys \
      -o results/wp.out -xi yes

   return

}



# ----------------------------------------------------------------------------------- #
# main
# ----------------------------------------------------------------------------------- #

main $@
