#
# TODO:
#

CC := mpic++ # Replace with g++ if not using mpi, also comment MPIFLAG
# CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR := test
BUILDDIR := build
TARGETDIR := bin
TARGET := bin/unittest

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -O3 -Wall -std=c++14 -DNDEBUG -flto
LIB :=
TESTOBJ := $(BUILDDIR)/catchtest.o
INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

unittest: $(BUILDDIR)/$(fname).o $(TESTOBJ)
	@echo " making $(fname)..."
	@echo " $(CC) $(BUILDDIR)/$(fname).o $(CFLAGS) $(INC) -o $(TARGETDIR)/$(fname) $(TESTOBJ) $(LIB)"; $(CC) $(BUILDDIR)/$(fname).o $(CFLAGS) $(INC) -o $(TARGETDIR)/$(fname) $(TESTOBJ) $(LIB)

example:
	@echo " making $(fname)..."
	@echo " $(CC) examples/$(fname).cpp $(CFLAGS) $(INC) -o $(TARGETDIR)/$(fname) $(LIB)"; $(CC) examples/$(fname).cpp $(CFLAGS) $(INC) -o $(TARGETDIR)/$(fname) $(LIB)

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Tests
# tester:
#   $(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester
#
#   # Spikes
#ticket:
#	  $(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean
