include ../../var.mk

LIB_FILE = libwf.a

all: $(LIB_FILE)

$(LIB_FILE): bohm-wavefunction-on-box-1d.o
	$(AR) rcs -o $@ $^

bohm-wavefunction-on-box-1d.o: ./bohm-wavefunction-on-box-1d.cc
	$(CXX) -o $@ $< -c -fPIC $(OPT_FLAGS) $(WARN_FLAGS) \
		$(TDSE_INC_FLAGS) $(TDSE_LIB_FLAGS)

clean:
	$(RM) *.o *.a

