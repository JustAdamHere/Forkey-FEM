BIN_DIR  ?= ./bin
LIB_DIR  ?= ./lib
SRC_DIR  ?= ./src
TEST_DIR ?= ./tests
TEST_FILE ?= test

include makefile.sources
OBJS := $(SRCS:%=$(BIN_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIR) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

FC     = gfortran
FFLAGS = -fcheck=all -fbacktrace -Wall -g

all:
	(make -f makefile $(OBJS))

lib:
	(make -f makefile all)
	(ar -r $(LIB_DIR)/libforkeyfem.a $(OBJS))

$(BIN_DIR)/%.out: $(TEST_DIR)/%.f90
	(make -f makefile all)
	#$(FC) $(FFLAGS) $< -I $(BIN_DIR) -o $@
	$(FC) $(FFLAGS) $(OBJS) $< -I $(BIN_DIR) -o $@

$(BIN_DIR)/%.f90.o: $(SRC_DIR)/%.f90
	mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -J $(BIN_DIR) -c $< -o $@

# $(INC_DIR)/%.mod: $(SRC_DIR)/%.f90
# 	mkdir -p $(dir $@)
# 	$(FC) $(FFLAGS) -c $< -o $@

clean ::
	$(RM) -r $(BIN_DIR)/*.mod   # Module files.
	$(RM) -r $(BIN_DIR)/*.out   # Test files.
	$(RM) -r $(BIN_DIR)/*.f90.o # Object files.
	$(RM) -r $(LIB_DIR)/*.a     # Library files.

.PHONY: clean