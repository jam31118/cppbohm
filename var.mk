# Requires the following environment variables:
# - `LAPACKHOME`
# - `GSLHOME`

OPT_FLAGS = -O3
WARN_FLAGS = -Wall

LAPACK_LIB_FLAGS = -L$(LAPACKHOME) -llapack -lrefblas -lgfortran

GSL_LIB_PATH = $(GSLHOME)/lib
GSL_LIB_FLAGS = -L$(GSL_LIB_PATH) -lgsl -lgslcblas -lm -Wl,-rpath,$(GSL_LIB_PATH)
GSL_INC_FLAGS = -I$(GSLHOME)/include

