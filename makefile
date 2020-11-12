BUILD_DIR ?= ./build
SRC_DIRS ?= ./src
TEST_DIRS ?= ./tests
TEST_FILE ?= test

TARGET_EXEC ?= $(TEST_FILE).out

SRCS := $(wildcard $(SRC_DIRS)/*.f90 -or -wholename $(TEST_DIRS)/$(TEST_FILE).f90)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

FCFLAGS := $(INC_FLAGS)
FC = gfortran

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	#$(MKDIR_P) $(dir $@)
	$(FC) $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.f90.o: %.f90
	$(MKDIR_P) $(dir $@)
	$(FC) $(FCFLAGS) -J $(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/%.mod: %.f90
	$(MKDIR_P) $(dir $@)
	$(FC) $(FCFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p