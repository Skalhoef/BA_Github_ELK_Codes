
MAKE = make
AR = ar


# Elk requires BLAS/LAPACK and fast Fourier transform (FFT) libraries.

# We recommend using optimised, mutithreaded BLAS/LAPACK libaries such as those
# provided by the Intel Math Kernel Library (MKL), OpenBLAS or BLIS.

# Elk can use either the Fastest Fourier Transform in the West (FFTW) or the
# MKL FFT. Both single and double precision versions of the FFT are needed.


#------------------------------------------------------------------------------#
#                                  Libraries                                   #
#------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------
SRC_MKL = mkl_stub.f90
# To enable MKL multithreaded parallelism, uncomment the following line and link
# with the MKL library.
#SRC_MKL =
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# If your compiler does *not* support the Message Passing Interface (MPI) then
# uncomment the line below.
#SRC_MPI = mpi_stub.f90
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Elk supports two FFT libraries:
#   Fastest Fourier Transform in the West (FFTW)
#   Intel Math Kernel Library (MKL) Fast Fourier Transform
#
# FFTW is enabled by default with:
SRC_FFT = zfftifc_fftw.f90 cfftifc_fftw.f90
#
# To enable MKL FFT instead, copy mkl_dfti.f90 to the elk/src directory and
# uncomment the line below.
#SRC_FFT = mkl_dfti.f90 zfftifc_mkl.f90 cfftifc_mkl.f90
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
SRC_LIBXC = libxcifc_stub.f90
# To enable Libxc first download and compile version 6.x of the library. Next
# copy the files libxcf90.a and libxc.a to the elk/src directory and uncomment
# the following lines.
#LIB_LIBXC = libxcf90.a libxc.a
#SRC_LIBXC = libxcf90.f90 libxcifc.f90
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
SRC_W90S = w90_stub.f90
# To enable the Wannier90 library copy libwannier.a to the elk/src directory and
# uncomment the following lines.
#SRC_W90S =
#LIB_W90 = libwannier.a
#-------------------------------------------------------------------------------



#------------------------------------------------------------------------------#
#                                  Compilers                                   #
#------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------
# Intel MPI Fortran version 18 and later.
#F90 = mpiifort
#F90_OPTS = -O3 -xHost -ipo -qopenmp -mkl=parallel
#F90_LIB = -liomp5 -lpthread -lm -ldl
#SRC_MKL =
#AR = xiar
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# GNU Fortran compiler with MKL.
#F90 = mpif90
#F90_OPTS = -Ofast -march=native -mtune=native -fopenmp -ffpe-summary=none -Wno-lto-type-mismatch
#F90_LIB = -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread
#SRC_MKL =
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# GNU Fortran compiler with BLIS, libflame and FFTW
#F90 = mpif90
#F90_OPTS = -Ofast -march=native -mtune=native -fopenmp -ffpe-summary=none -Wno-lto-type-mismatch
#F90_LIB = -lblis -lflame -lfftw3 -lfftw3f
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# GNU Fortran compiler with OpenBLAS, LAPACK and FFTW
F90 = mpif90
F90_OPTS = -Ofast -march=native -mtune=native -fopenmp -ffpe-summary=none -Wno-lto-type-mismatch
F90_LIB = -lopenblas -llapack -lfftw3 -lfftw3f
SRC_OBLAS =
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Intel Fortran with debugging and profiling options.
#F90 = mpiifort
#F90_OPTS = -O3 -qopenmp -mkl=parallel -pg -xHost -init=snan,arrays -traceback
#F90_LIB = -liomp5 -lpthread -lm -ldl
#SRC_MKL =
#AR = ar
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# GNU Fortran debugging and profiling options.
#F90 = gfortran
#F90_OPTS = -O3 -Wall -Wunused -Warray-temporaries -Waliasing -Wampersand -Wconversion -Winteger-division -Wintrinsics-std -Wsurprising -Wtabs -Wunused-parameter -Wfunction-elimination -Wrealloc-lhs-all -Wcompare-reals -fopenmp -fbounds-check -pg
#F90_LIB = -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread
#SRC_MKL =
#SRC_MPI = mpi_stub.f90
#-------------------------------------------------------------------------------

