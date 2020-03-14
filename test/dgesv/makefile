include ../../var.mk

bin = dgesv

BOHM_HOME = ../..
BOHM_LIB_DIR = $(BOHM_HOME)/lib
BOHM_LIB_FLAGS = '-Wl,-rpath,$$ORIGIN/$(BOHM_LIB_DIR)' -L$(BOHM_LIB_DIR) -lbohm

all: $(bin)

$(bin): libbohm.so
	$(CXX) -o $@ $@.cc $(OPT_FLAGS) $(WARN_FLAGS) $(LAPACK_LIB_FLAGS) $(BOHM_LIB_FLAGS)

libbohm.so:
	$(MAKE) -C $(BOHM_HOME)

clean:
	$(RM) $(bin)

