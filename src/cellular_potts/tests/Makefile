# Default target, for when you just run make
.PHONY: test
test: run_all_tests


# Get includes and libraries for Catch2
# We skip this when doing make clean, because we don't need the information and
# Catch2 may not be available, which would cause this to error out.
ifneq "$(filter $(MAKECMDGOALS),clean)" "clean"
    PCPATH := $(PKG_CONFIG_PATH):../../../lib/Catch2/catch2/share/pkgconfig
    CATCH2_INCLUDES := $(shell PKG_CONFIG_PATH=$(PCPATH) pkg-config --cflags catch2-with-main)
    CATCH2_LIBS := $(shell PKG_CONFIG_PATH=$(PCPATH) pkg-config --libs catch2-with-main)

    CXXFLAGS := $(CATCH2_INCLUDES) $(CXXFLAGS) -g
    CXXFLAGS += -I. -I.. -I../.. -I../../graphics -I../../models
    CXXFLAGS += -I../../parameters -I../../plotting -I../../reaction_diffusion
    CXXFLAGS += -I../../util -I../../xpm -I../../compute
    CXXFLAGS += -I../../spatial -I../../parameters
    CXXFLAGS += -I../../../lib/MultiCellDS/v1.0/v1.0.0/libMCDS/mcds_api/
    CXXFLAGS += -I../../../lib/MultiCellDS/v1.0/v1.0.0/libMCDS/xsde/libxsde
    LDFLAGS := $(CATCH2_LIBS) $(LDFLAGS)

    CATCH2_INCLUDE_DIR := ../../../lib/Catch2/catch2/include
endif

# Find tests by name, then remove the .cpp extension
TESTS := $(patsubst %.cpp, %, $(wildcard test_*.cpp))
TEST_EXECUTABLES := $(patsubst %,build/%, $(TESTS))


# Define targets that run tests
.PHONY: run_%
run_%: build/%
	./$^

# List all the run-a-test targets and create a target depending on them all.
# We include the test executables explicitly here, or Make will consider them
# intermediate targets and remove them at the end of the run!
RUN_TARGETS := $(patsubst %,run_%,$(TESTS))

.PHONY: run_all_tests
run_all_tests: $(TEST_EXECUTABLES) $(RUN_TARGETS)


# Find dependencies for the tests, so that they get rebuilt if you change any
# headers they include. Note that dependencies on source files still need to
# be specified by hand, and that if you change which headers are included by
# a header, you need to make clean and rebuild from scratch.
#
# The C++ compiler, when given the -MM option and a file, will scan all the
# included headers and produce output in Make format specifying the
# dependencies. We save that to a file with a .d extension and the same name
# as the test. We mark the Catch2 include directory as as system directory so
# that -MM will not include any Catch2 headers in the output.
build/test_%.d: test_%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -isystem $(CATCH2_INCLUDE_DIR) -E -MM -MT $(@:.d=) -MF $@ $<

# If you try to include a file that does not exist, Make will try to build it,
# in this case using the rule above. We don't include dependencies if we're
# running "make clean", because that would build them and we're actually trying
# to clean up.
ifneq "$(filter $(MAKECMDGOALS),clean)" "clean"
    DEPS := $(TESTS:%=build/%.d)
    include $(DEPS)
endif

build/test_%: test_%.cpp
	$(CXX) -o $@ $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS)


clean:
	rm -f $(TEST_EXECUTABLES) build/*.d
