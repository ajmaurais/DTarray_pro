
#
# C++ Compiler
CXX := g++
#
#
# Flags
#
#   Compiler
CXXFLAGS += -c -g -Wall -std=c++11
#
#   Linker
LDFLAGS += -g
#
#   Library
LIBFLAGS :=
#
#   Include
INCLUDEFLAGS :=
#
#
# Program name
EXE := DTarray_pro
#
#
# Directories
#
#   Headers
HEADERDIR := include
#
#   Sources
SRCDIR := src
#
#   Objects
OBJDIR := obj
#
#   Binary
BINDIR := bin
#
#
#   Install dirrectory
INSTALL_DIR := /usr/local/bin/
#
#
################################################################################

HEADERS := $(wildcard $(HEADERDIR)/*.h)
SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(subst $(SRCDIR)/,$(OBJDIR)/,$(SRCS:.cpp=.o))

CXXFLAGS += $(INCLUDEFLAGS) -I$(HEADERDIR)
LDFLAGS += $(LIBFLAGS)

.PHONY: all clean distclean install

all: $(BINDIR)/$(EXE) helpFile.pdf

$(BINDIR)/$(EXE): $(OBJS)
	mkdir -p $(BINDIR)
	$(CXX) $(LDFLAGS) $? -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

helpFile.pdf : db/helpFile.man
	bash updateMan.sh

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(EXE)
	rm -f helpFile.pdf

install: $(BINDIR)/$(EXE)
	cp $(BINDIR)/$(EXE) $(INSTALL_DIR)$(EXE)

distclean: clean
