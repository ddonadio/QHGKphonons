FC=mpif90
FFLAGS=" -fopenmp -m64"
LAPACK="-lmkl_lapack95_lp64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread"

#Create bin folder
mkdir -p ../bin

${FC} ${FFLAGS} -o ../bin/relaxmpi.x relaxmpi.f90 ${LAPACK}
${FC} ${FFLAGS} -o ../bin/RdiagDynMat.x diagDynMat.f90 ${LAPACK}
${FC} ${FFLAGS} -o ../bin/QHGKd.x QHGKd.f90 ${LAPACK}
