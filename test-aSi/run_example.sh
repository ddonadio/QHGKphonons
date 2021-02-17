#!/bin/bash

DIAG_EXE="../../bin/RdiagDynMat.x"
RMPI="../../bin/relaxmpi.x"
QHGK="../../bin/QHGKd.x"
THREADS=4
PROCS=1

export MKL_NUM_THREADS=$THREADS
natom=512
temperature=300
delta=0.5
nspecies=1

gunzip THIRD.gz Dyn.form.gz

mkdir -p results
cp CONFIG CONTROL FIELD Dyn.form REFERENCE THIRD results/
cd results

$DIAG_EXE ${natom} -sym
$RMPI ${natom} ${temperature} ${delta}
tail -1533 Decay.000 > DECAY
$QHGK ${natom} ${nspecies} -bX -bY -bZ
