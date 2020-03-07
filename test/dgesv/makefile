include ../var.mk

all: dgesv

dgesv:
	$(CXX) -o $@ $@.cc -I$(BOHMHOME)/include $(LAPACK-LIB-FLAGS)

clean:
	$(RM) dgesv

