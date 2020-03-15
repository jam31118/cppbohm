# Requires the following environment variables:
# - `LAPACKHOME`
# - `GSLHOME`

OPT_FLAGS = -O3
WARN_FLAGS = -Wall

LAPACK_LIB_RPATH_FLAGS = -Wl,-rpath,$(LAPACKHOME)
LAPACK_LIB_FLAGS = -L$(LAPACKHOME) -llapack -lrefblas -lgfortran

GSL_LIB_PATH = $(GSLHOME)/lib
GSL_LIB_FLAGS = -L$(GSL_LIB_PATH) -lgsl -lgslcblas -lm -Wl,-rpath,$(GSL_LIB_PATH)
GSL_INC_FLAGS = -I$(GSLHOME)/include

TDSE_INC_FLAGS = -I$(TDSEHOME)/include
TDSE_LIB_PATH = $(TDSEHOME)/lib
TDSE_LIB_FLAGS = -Wl,-rpath,$(TDSE_LIB_PATH) -L$(TDSE_LIB_PATH) -ltdse

PARAM_INC_FLAGS = -I$(PARAMHOME)/include
PARAM_LIB_PATH = $(PARAMHOME)/lib
PARAM_LIB_FLAGS = -Wl,-rpath,$(PARAM_LIB_PATH) -L$(PARAM_LIB_PATH) -lparam

