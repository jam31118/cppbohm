include ../../var.mk

bin = rosenbrock

all: $(bin)

$(bin): $(bin).cc
	$(CXX) -o $@ $< $(OPT_FLAGS) $(WARN_FLAGS) $(GSL_INC_FLAGS) $(GSL_LIB_FLAGS)

clean:
	$(RM) $(bin)


