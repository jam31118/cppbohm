include ../../var.mk

LIB_FILE = libpropagator.a

all: $(LIB_FILE)

$(LIB_FILE): bohm-propagator-on-box-1d.o bohm-propagator-on-box-1d-with-field.o
	$(AR) rcs -o $@ $^

bohm-propagator-on-box-1d.o: bohm-propagator-on-box-1d.cc
	$(CXX) -o $@ $< -c -fPIC $(OPT_FLAGS) $(WARN_FLAGS) \
		$(TDSE_INC_FLAGS) \
		$(GSL_INC_FLAGS) \
		$(DEFINE_FLAGS)

bohm-propagator-on-box-1d-with-field.o: bohm-propagator-on-box-1d-with-field.cc
	$(CXX) -o $@ $< -c -fPIC $(OPT_FLAGS) $(WARN_FLAGS) \
		$(TDSE_INC_FLAGS) \
		$(GSL_INC_FLAGS)

clean:
	$(RM) *.o *.a

