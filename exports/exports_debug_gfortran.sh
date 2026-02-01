FPM_FC="gfortran"
export FPM_FC

FPM_FFLAGS_DEBUG="-DDEBUG -g -fmax-errors=1 -ffixed-line-length-none -Wall -Wimplicit-interface -Wall -Wno-integer-division -fbounds-check"
FPM_FFLAGS="$FPM_FFLAGS_DEBUG"
export FPM_FFLAGS
