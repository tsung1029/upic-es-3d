#!/bin/ksh
#
# Mnew_pbeps32.ksh -- atlas parallel run script
#
#MSUB -S /bin/ksh
#MSUB -r n                          ##no restart
#MSUB -o
#MSUB -j oe                         ##combine error and output logs

set -xv
umask 027

##############################################################################
#test variables: needs HOME, LOGNAME, RUNDIR, PFS_TOPDIR from the submit line

if [ "$RUNDIR"= "" ]; then
   echo "cannot run with RUNDIR"
   exit 1
fi

cd ${RUNDIR}

##############################################################################
# perform the run

if [ "$SLURM_NPROCS" = "" ]; then
   typeset -x -i SLURM_NPROCS=8*SLURM_NNODES
fi

# in the below line
#     -l     label output by prepending process rank number
#     -i0    connect stdin only to rank 0
#     -c1    one cpu per task (ie, single threading)
#

srun -l -i0 -c1 -n $SLURM_NPROCS new_pbeps32.out

##############################################################################
# All done

exit

#MSUB -m

