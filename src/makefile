include ../var.mk

LIB_DIR = ../lib
LIB_FILE = libbohm.so

all: install

install: $(LIB_FILE)
	if [ ! -d "$$LIB_DIR" ]; then mkdir $(LIB_DIR); fi; mv $(LIB_FILE) $(LIB_DIR)

$(LIB_FILE): lapack.o
	$(CXX) -o $@ $^ -shared

lapack.o: lapack.cc
	$(CXX) -o $@ $< -c -fPIC $(OPT_FLAGS) $(WARN_FLAGS)

clean:
	$(RM) *.o *.so $(LIB_DIR)/$(LIB_FILE)
