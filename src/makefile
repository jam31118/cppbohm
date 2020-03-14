include ../var.mk

LIB_DIR = ../lib
LIB_FILE = libbohm.so

WF_LIB = libwf.a
WF_LIB_DIR = ./wf
PROPAGATOR_LIB = libpropagator.a
PROPAGATOR_LIB_DIR = ./propagator

all: install

install: $(LIB_FILE)
	if [ ! -d "$(LIB_DIR)" ]; then mkdir $(LIB_DIR); fi; mv $^ $(LIB_DIR)

$(LIB_FILE): lapack.o $(WF_LIB) $(PROPAGATOR_LIB)
	$(CXX) -o $@ -shared *.o \
		-Wl,--whole-archive \
		$(WF_LIB_DIR)/$(WF_LIB) $(PROPAGATOR_LIB_DIR)/$(PROPAGATOR_LIB) \
		-Wl,--no-whole-archive

lapack.o: lapack.cc
	$(CXX) -o $@ $< -c -fPIC $(OPT_FLAGS) $(WARN_FLAGS)

$(WF_LIB):
	$(MAKE) -C $(WF_LIB_DIR)

$(PROPAGATOR_LIB):
	$(MAKE) -C $(PROPAGATOR_LIB_DIR)

clean:
	$(RM) *.o *.so $(LIB_DIR)/$(LIB_FILE)
	$(MAKE) clean -C $(WF_LIB_DIR)
	$(MAKE) clean -C $(PROPAGATOR_LIB_DIR)

