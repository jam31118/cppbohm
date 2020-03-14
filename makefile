LIB_DIR = ./lib
SRC_DIR = ./src
TEST_DIR = ./test

all: install

install:
	$(MAKE) -C $(SRC_DIR)

clean:
	$(MAKE) clean -C $(SRC_DIR)
	$(MAKE) clean -C $(TEST_DIR)
	if [ -d "$(LIB_DIR)" ]; then rmdir $(LIB_DIR); fi
