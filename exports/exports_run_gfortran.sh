FPM_FC="gfortran"
export FPM_FC

FPM_FFLAGS_RUN="-O3 -march=native -fmax-errors=1 -ffixed-line-length-none -fno-range-check"
FPM_FFLAGS="$FPM_FFLAGS_RUN"
export FPM_FFLAGS

OMP_NUM_THREADS=8
export OMP_NUM_THREADS
