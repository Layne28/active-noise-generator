#Some organization borrowed from: https://github.com/tscott8706/cpp-csv-col-replacer/blob/master

# Build executable with:
# % make
# Delete object files and executable with:
# % make clean
# Rebuild all objects and executable with:
# % make -B

SRC_DIR := src
OBJ_DIR := build
BIN_DIR := bin
LIB_DIR := lib
LIB_OBJ_DIR := lib/build
TEST_SRC_DIR := test/src
TEST_OBJ_DIR := test/build
ANALYSIS_SRC_DIR := analysis/src
ANALYSIS_OBJ_DIR := analysis/build
LIB := $(LIB_DIR)/libangen.so

EXECUTABLE := active_noise_generator
TEST_EXECUTABLE := test_active_noise_generator
ANALYSIS_EXECUTABLE := analyze
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
TEST_SOURCES := $(wildcard $(TEST_SRC_DIR)/*.cpp)
ANALYSIS_SOURCES := $(wildcard $(ANALYSIS_SRC_DIR)/*.cpp)
HEADERS := $(wildcard $(SRC_DIR)/*.hpp)
TEST_HEADERS := $(wildcard $(TEST_SRC_DIR)/*.hpp)
ANALYSIS_HEADERS := $(wildcard $(ANALYSIS_SRC_DIR)/*.hpp)

CXX := g++

SHELL = /bin/sh

# Flags to pass to the compiler; per the reccomendations of the GNU Scientific Library
CXXFLAGS:= -std=c++17 -Wextra -pedantic -Wall -W -Wmissing-declarations -Wuninitialized -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -fshort-enums -fno-common -m64 -fopenmp -I$(HOME)/.local/include

# Compiler flags controling optimization levels. Use -O3 for full optimization,
# but make sure your results are consistent
# -g includes debugging information. You can also add -pg here for profiling 
PROFILE=-pg
OPTFLAGS:=$(PROFILE) -O2 #Might try changing to O3 to increase speed

# Flags to pass to the linker; -lm links in the standard c math library
#LDFLAGS:= -fopenmp -lfftw3 -lm -lgsl -lgslcblas -llapack -lblas -larmadillo -lstdc++fs $(PROFILE) -L$(HOME)/.local/lib 
LDFLAGS:= -fopenmp -lfftw3 -lm -lgsl -lgslcblas -lopenblas -larmadillo -lstdc++fs $(PROFILE) -L$(HOME)/.local/lib 

# Variable to compose names of object files from the names of sources
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
OBJECTS_NO_MAIN = $(filter-out $(OBJ_DIR)/main.o,$(OBJECTS))
OBJECTS_LIB := $(patsubst $(SRC_DIR)/%.cpp,$(LIB_OBJ_DIR)/%.o,$(SOURCES_NO_MAIN))

#When compiling tests, include all objects in actual program except for main
#(there's a main function in the test folder)
TEST_OBJECTS := $(patsubst $(TEST_SRC_DIR)/%.cpp,$(TEST_OBJ_DIR)/%.o,$(TEST_SOURCES))
TEST_OBJECTS += $(OBJECTS_NO_MAIN)

ANALYSIS_OBJECTS := $(patsubst $(ANALYSIS_SRC_DIR)/%.cpp,$(ANALYSIS_OBJ_DIR)/%.o,$(ANALYSIS_SOURCES))

# Default target depends on sources and headers to detect changes
all: $(SOURCES) $(HEADERS)  $(BIN_DIR)/$(EXECUTABLE)
analysis: $(ANALYSIS_SOURCES) $(ANALYSIS_HEADERS) $(BIN_DIR)/$(ANALYSIS_EXECUTABLE)
lib: $(LIB)
	install $(LIB) $(HOME)/.local/lib/
	mkdir -p $(HOME)/.local/include/libangen
	install $(HEADERS) $(HOME)/.local/include/libangen/
install: 
	install bin/* $(HOME)/.local/bin/
test: $(TEST_SOURCES) $(TEST_HEADERS) $(BIN_DIR)/$(TEST_EXECUTABLE)

# Rule to compile a source file to object code
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $(OPTFLAGS) $< -o $@
$(LIB_OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -fPIC -c $(CXXFLAGS) $(OPTFLAGS) $< -o $@
$(TEST_OBJ_DIR)/%.o : $(TEST_SRC_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@
$(ANALYSIS_OBJ_DIR)/%.o : $(ANALYSIS_SRC_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

# Build the executable by linking all objects
$(BIN_DIR)/$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@
$(LIB): $(OBJECTS_LIB)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -o $@ $(OBJECTS_LIB)
$(BIN_DIR)/$(TEST_EXECUTABLE): $(TEST_OBJECTS)
	$(CXX) $(TEST_OBJECTS) $(LDFLAGS) -o $@
$(BIN_DIR)/$(ANALYSIS_EXECUTABLE): $(ANALYSIS_OBJECTS)
	$(CXX) $(ANALYSIS_OBJECTS) $(LDFLAGS) -o $@

# clean up so we can start over (removes executable!)
clean:
	rm -f $(OBJ_DIR)/*.o $(LIB_OBJ_DIR)/*.o $(TEST_OBJ_DIR)/*.o $(ANALYSIS_OBJ_DIR)/*.o $(BIN_DIR)/$(EXECUTABLE) $(BIN_DIR)/$(TEST_EXECUTABLE) $(BIN_DIR)/$(ANALYSIS_EXECUTABLE)
