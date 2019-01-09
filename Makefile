CXX ?= g++

# path #
SRC_PATH = .
BUILD_PATH = build
BIN_PATH = $(BUILD_PATH)/bin

# executable # 
BIN_NAME = ppm1d

# extensions #
SRC_EXT = cpp

# code lists #
# Find all source files in the source directory, sorted by
# most recently modified
SOURCES = $(shell find $(SRC_PATH) -name '*.$(SRC_EXT)' | sort -k 1nr | cut -f2-)
# Set the object file names, with the source directory stripped
# from the path, and the build path prepended in its place
OBJECTS = $(SOURCES:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)
# Set the dependency files that will be used to add header dependencies
DEPS = $(OBJECTS:.o=.d)

# compile variables #
LOOP_PAR_THREADS = 1
# flags common #
FLTO_FLAGS 		= #-flto -fuse-linker-plugin
OPENMP_FLAGS 	= #-fopenmp
ILS_FLAGS 		= -ftree-vectorize -floop-interchange -floop-strip-mine \
								-floop-block -ftree-loop-distribution \
								-ftree-loop-distribute-patterns
STD_FLAGS 		= -g -std=c++17 -fdiagnostics-color $(FLTO_FLAGS) $(OPENMP_FLAGS)

# sanatize flags#
SANATIZE_FLAGS = -fsanitize=address -fsanitize=alignment -fstack-check -fsanitize=undefined -fsanitize=float-divide-by-zero
# debug flags #
DBG_FLAGS = -g3 -Wall -Wextra -Og -ggdb -fmax-errors=5 -fno-omit-frame-pointer -fno-inline $(SANATIZE_FLAGS)
# opt flags #
OPT_FLAGS = -O2 -march=native -Winline
#XPT_FLAGS = -O3 -march=native
XPT_FLAGS = -O3 -march=native -mtune=native $(ILS_FLAGS) -Winline
#XPT_FLAGS = -g -Ofast -ffast-math -march=native

#COMPILE_FLAGS_PROFILE = $(STD_FLAGS) -fprofile-use
#COMPILE_FLAGS_PROFILE = $(STD_FLAGS) -fprofile-generate

COMPILE_FLAGS_FIX = $(STD_FLAGS) $(DBG_FLAGS)
COMPILE_FLAGS_OPT = $(STD_FLAGS) $(OPT_FLAGS)
COMPILE_FLAGS_FAST = $(STD_FLAGS) $(XPT_FLAGS)


INCLUDES = -I./include/ -I/usr/include -I/$(HOME)/.local/include

# linker
LINK = $(CXX)
# linker flags
LD_FLAGS = 
# Space-separated libraries
LIBS = -lm 

.PHONY: default_target
default_target: debug

.PHONY: debug
debug: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS_FIX)
debug: dirs
	@$(MAKE) all

.PHONY: opt
opt: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS_OPT)
opt: dirs
	@$(MAKE) all

.PHONY: fast
fast: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS_FAST)
fast: dirs
	@$(MAKE) all

.PHONY: profile
profile: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS_PROFILE)
profile: dirs
	@$(MAKE) all

.PHONY: dirs
dirs:
	@echo "Creating directories"
	@mkdir -p $(dir $(OBJECTS))
	@mkdir -p $(BIN_PATH)

.PHONY: clean
clean:
	@echo "Deleting $(BIN_NAME) symlink"
	@$(RM) $(BIN_NAME)
	@echo "Deleting directories"
	@$(RM) -r $(BUILD_PATH)
	@$(RM) -r $(BIN_PATH)

# checks the executable and symlinks to the output
.PHONY: all
all: $(BIN_PATH)/$(BIN_NAME)
	@echo "Making symlink: $(BIN_NAME) -> $<"
	@$(RM) $(BIN_NAME)
	@ln -s $(BIN_PATH)/$(BIN_NAME) $(BIN_NAME)

# Creation of the executable
$(BIN_PATH)/$(BIN_NAME): $(OBJECTS)
	@echo "Linking: $@"
	$(LINK) $(OBJECTS) -o $@ $(LD_FLAGS) $(CXXFLAGS) $(LIBS)

# Add dependency files, if they exist
-include $(DEPS)

# Source file rules
# After the first compilation they will be joined with the rules from the
# dependency files to provide header dependencies
$(BUILD_PATH)/%.o: $(SRC_PATH)/%.$(SRC_EXT)
	@echo "Compiling: $< -> $@"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@
