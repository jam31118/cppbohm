include ../../../var.mk

BOHM_HOME = $(CURDIR)/../../..
BOHM_LIB_DIR = $(BOHM_HOME)/lib
BOHM_LIB_FLAGS = -Wl,-rpath,$(BOHM_LIB_DIR) -L$(BOHM_LIB_DIR) -lbohm

bin = box-1d
bin_field = box-1d-with-field

all: $(bin) $(bin_field)

$(bin): $(bin).cc libbohm.so
	$(CXX) -o $@$(SUFFIX) $< $(OPT_FLAGS) $(WARN_FLAGS) \
		$(TDSE_INC_FLAGS) $(PARAM_INC_FLAGS) \
		$(LAPACK_LIB_RPATH_FLAGS) $(LAPACK_LIB_FLAGS) \
		$(BOHM_LIB_FLAGS) $(TDSE_LIB_FLAGS) $(PARAM_LIB_FLAGS) \
		$(GSL_INC_FLAGS) $(GSL_LIB_FLAGS) \
		$(DEF_FLAGS)

$(bin_field):
	$(MAKE) $(bin) DEF_FLAGS=-DFIELD SUFFIX=-with-field

libbohm.so:
	$(MAKE) -C $(BOHM_HOME)/src

clean:
	$(RM) $(bin) $(bin_field)
