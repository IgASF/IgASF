#DISABLE STANDARD RULES
MAKEFLAGS += --no-builtin-rules
.SUFFIXES:

#TYPE OF BUILD AND PERSONALIZED FLAGS
TYPE ?=debug

#PLACING OF OUTPUT
BLDDIR:=build/$(TYPE)
OBJDIR:=$(BLDDIR)/obj
DEPDIR:=$(BLDDIR)/dep
BINDIR:=$(BLDDIR)/bin
$(shell mkdir -p $(DEPDIR) >/dev/null)
$(shell mkdir -p $(OBJDIR) >/dev/null)
$(shell mkdir -p $(BINDIR) >/dev/null)

#SRC DIRECTORIES
LIB_DIRS    := src/algebra/ src/assembling/ src/bases/ src/quadrature/ src/tools/ src/maps/ src/solve/
BIN_DIRS    := src/bin/

#COMPILING FLAGS
#$(foreach dir, $(LIB_DIRS), -I $(dir))
INCLUDE_FLAGS:= -I extern/eigen -I extern/json -I extern/unittests  -I src
DEPFLAGS= -MD -MP -MF $(DEPDIR)/$(1).dep -MT $$@
export INCLUDE_FLAGS
export DEPFLAGS

# Create those files if they do not yet exist
$(shell if [[ ! -e config/debug ]]; then cp config/debug_template config/debug; fi)
$(shell if [[ ! -e config/release ]]; then cp config/release_template config/release; fi)

# Give some defaults
CXXFLAGS:= -O0 -g -fPIC $(INCLUDE_FLAGS) -Wextra -Wall -Werror -std=c++1z -ftemplate-backtrace-limit=0 -pthread
ECHO_FLAG:= -e
LIB_EXT:= so

-include config/$(TYPE)

#FILES BY TYPE
LIB_SRC     :=$(wildcard $(addsuffix *.cpp, $(LIB_DIRS)))
LIB_OBJ     :=$(patsubst %.cpp, $(OBJDIR)/%.o, $(notdir $(LIB_SRC)))
LIB         :=$(BLDDIR)/libassemble.$(LIB_EXT)

BIN_SRC     :=$(wildcard $(addsuffix *.cpp, $(BIN_DIRS)))
BIN         :=$(patsubst %.cpp, $(BINDIR)/%, $(notdir $(BIN_SRC)))
EXE_LN      :=$(patsubst %.cpp, %-$(TYPE), $(notdir $(BIN_SRC)))

SCR_SRC     :=$(wildcard $(addsuffix *.sh, $(BIN_DIRS)))
SCR         :=$(patsubst %.sh, $(BINDIR)/%.sh, $(notdir $(SCR_SRC)))
SCR_LN      :=$(patsubst %.sh, %-$(TYPE), $(notdir $(SCR_SRC)))

# unittests and test-cases
UTS_SRC     :=$(wildcard $(addsuffix *.uts,  $(LIB_DIRS)))
UTS_OBJ     :=$(patsubst %.uts, $(OBJDIR)/%.uts.o, $(notdir $(UTS_SRC)))

UTS_MAIN    :=$(BIN_DIRS)unittest.cc
UTS         :=$(BINDIR)/unittest
UTS_LN      :=unittest-$(TYPE)


#TARGETS
.PHONY: help
help:
	@echo "                    _____                ________________                   "
	@echo "                   |_   _|        /\    / ______________/                   "
	@echo "                     | |  __ _   /  \  | (___ |  ____/                      "
	@echo "                     | | / _' | / /\ \  \___ \| |__                         "
	@echo "                    _| || (_| |/ ____ \     ) |  __|                        "
	@echo "                   |_____\__, /_/   _\_\____) | |                           "
	@echo "                          __/ |   /__________/|_|                           "
	@echo "                         |___/                                              "
	@echo "                                                                            "
	@echo " -------------------------------------------------------------------------- "
	@echo "                Isogeometric-Assembling-by-Sum-Factorization                "
	@echo " -------------------------------------------------------------------------- "
	@echo "                                                                            "
	@echo " https://github.com/IgASF/IgASF                                             "
	@echo "                                                                            "
	@echo " Type                                                                       "
	@echo "                                                                            "
	@echo "   make all               ... to make executables of this library (debug)   "
	@echo "   make all-release       ... to make executables of this library (release) "
	@echo "   make so                ... to make this library as so-file (debug)       "
	@echo "   make so-release        ... to make this library as so-file (release)     "
	@echo "   make unittest          ... to make and run unittests (debug)             "
	@echo "   make unittest-release  ... to make and run unittests (release)           "
	@echo "   make clean             ... to clean the build environment (debug)        "
	@echo "   make clean-release     ... to clean the build environment (release)      "
	@echo "   make echoall           ... to get the contents of the given variables    "
	@echo "                                                                            "
	@echo " Use the files                                                              "
	@echo "                                                                            "
	@echo "   config/debug   and   config/release                                      "
	@echo "                                                                            "
	@echo " to specify the corresponding configuration. If those files are deleted,    "
	@echo " they will be restored based on the corresponding templates.                "
	@echo "                                                                            "
	@echo " Invoke after 'make all'                                                    "
	@echo "                                                                            "
	@echo "   ./generateTest-debug   ... to create a configuration file for assembling "
	@echo "                              the problem of interst                        "
	@echo "                                                                            "
	@echo "   ./runTest-debug test   ... to assemble a problem specified by such a     "
	@echo "                              configuration file test, like with            "
	@echo "                                                                            "
	@echo "   ./runTest-debug tests/QuarterAnnulus_stiff_128/deg03 -o mat              "
	@echo "                                                                            "
	@echo "                              to assemble and store result in file mat using"
	@echo "                              a binary format                               "
	@echo "                                                                            "
	@echo "   ./echoMatrix-debug mat ... to echo the matrix data from file mat         "
	@echo "                                                                            "
	@echo "   ./compareMatrices-debug mat1 mat2                                        "
	@echo "                          ... to compare such matrix files                  "
	@echo "                                                                            "
	@echo "   ./runTestsAndCompare-debug test[s]                                       "
	@echo "                          ... to assemble the problem with all methods and  "
	@echo "                              compare the results and write results to      "
	@echo "                              log.txt and compare.txt, like with            "
	@echo "                                                                            "
	@echo "   ./runTestsAndCompare-debug tests/QuarterAnnulus_stiff_128/*              "
	@echo "                                                                            "
	

.PHONY: all so
all: $(EXE_LN) $(SCR_LN)
so: $(LIB)

.PHONY: echoall
echoall:
	@echo LIB_DIRS $(LIB_DIRS)
	@echo LIB_SRC $(LIB_SRC)
	@echo LIB_OBJ $(LIB_OBJ)
	@echo BIN_DIRS $(BIN_DIRS)
	@echo BIN_SRC $(BIN_SRC)
	@echo BIN $(BIN)
	@echo INCLUDE_FLAGS $(INCLUDE_FLAGS)
	@echo CXXFLAGS $(CXXFLAGS)
	@echo LDFLAGS $(CXXFLAGS)
	@echo DEPFLAGS $(DEPFLAGS)
	@echo UTS_MAIN $(UTS_MAIN)
	@echo UTS_SRC $(UTS_SRC)
	@echo SCR_SRC $(SCR_SRC)
	@echo SCR $(SCR)
	@echo SCR_LN $(SCR_LN)

#DEPENDENCIES
$(DEPDIR)/%.dep: ;
.PRECIOUS: $(OBJDIR)/%.o $(DEPDIR)/%.dep $(BLDDIR)/$(LIB_GSM)
.SECONDARY:

# functions

define make-object-rule
$(OBJDIR)/$(1).o: $(2) $(DEPDIR)/$(1).dep
	@echo $(ECHO_FLAG) "\033[0;31mmake $$@\033[0m"
	$(CXX) $(DEPFLAGS)  -c -o $$@ -x c++ $$<  $(CXXFLAGS)
endef

define make-bin-rule
$(BINDIR)/$(1): $(2) $(DEPDIR)/$(1).dep $(LIB_OBJ)
	@echo $(ECHO_FLAG) "\033[0;31mmake $$@\033[0m"
	$(CXX)  $(DEPFLAGS) -L $(BLDDIR) -o $$@ $$< $(CXXFLAGS) $(LIB_OBJ)
endef

define make-ln-rule
$(1)-$(TYPE): $(BINDIR)/$(1)
	@echo $(ECHO_FLAG) "\033[0;31mmake $$@\033[0m"
	ln -fs $$< $(1)-$(TYPE)
endef

define make-scr-rule
$(BINDIR)/$(1).sh: $(2)
	@echo $(ECHO_FLAG) "\033[0;31mmake $$@\033[0m"
	sed 's/TYPE/$(TYPE)/g' $$< > $$@
	chmod +x $$@
endef

define make-scr-ln-rule
$(1)-$(TYPE): $(BINDIR)/$(1).sh
	@echo $(ECHO_FLAG) "\033[0;31mmake $$@\033[0m"
	ln -fs $$< $(1)-$(TYPE)
endef

# library
$(foreach source, $(LIB_SRC),     $(eval $(call make-object-rule,$(basename $(notdir $(source))),$(source))))

$(LIB) : $(LIB_OBJ) $(DEPDIR)/%.dep
	@echo $(ECHO_FLAG) "\033[0;31mmake $@\033[0m"
	$(CXX) $(LDFLAGS) $(LIB_OBJ) --shared -o $@

# bins

$(foreach src, $(BIN_SRC),     $(eval $(call make-bin-rule,$(basename $(notdir $(src))),$(src))))
$(foreach src, $(BIN_SRC),     $(eval $(call make-ln-rule,$(basename $(notdir $(src))),$(src))))

$(foreach src, $(SCR_SRC),     $(eval $(call make-scr-rule,$(basename $(notdir $(src))),$(src))))
$(foreach src, $(SCR_SRC),     $(eval $(call make-scr-ln-rule,$(basename $(notdir $(src))),$(src))))

# unittests and checking
$(foreach src, $(UTS_SRC),  $(eval $(call make-object-rule,$(notdir $(src)),$(src))))

$(UTS) : $(UTS_OBJ) $(UTS_MAIN) $(LIB_OBJ)
	@echo $(ECHO_FLAG) "\033[0;31mmake $@\033[0m"
	$(CXX) $(CXXFLAGS) $(UTS_MAIN)  $(UTS_OBJ) -L $(BLDDIR) -o $(UTS) $(LIB_OBJ)

$(UTS_LN): $(UTS)
	@echo $(ECHO_FLAG) "\033[0;31mmake $@\033[0m"
	ln -fs $(UTS) $(UTS_LN)

.PHONY: unittest
unittest: $(UTS_LN)
	@echo $(ECHO_FLAG) "\033[0;31mmake $@\033[0m"
	./$(UTS_LN)

# cleaning rules

.PHONY: clean
clean:
	@echo $(ECHO_FLAG) "\033[0;31mmake $@\033[0m"
	-rm -rf build/$(TYPE) $(BIN) $(EXE_LN) $(SCR_LN) $(UTS_LN) || true
	-rm -f $(patsubst build/$(TYPE)/bin/%, %-$(TYPE), $(BIN) ) || true
	-rmdir build || true

# convenience

ifneq ($(TYPE), debug)
.PHONY: %-debug
%-debug:
	make TYPE=debug $(patsubst %-debug, %, $@)
endif

ifneq ($(TYPE), release)
.PHONY: %-release
%-release:
	make TYPE=release $(patsubst %-release, %, $@)
endif

-include $(DEPDIR)/*.dep
