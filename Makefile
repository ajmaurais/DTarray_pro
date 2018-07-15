
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
EXE := DTarray
#
#
# Directories
#
#   Headers
HEADERDIR := include
#
#   .git
GITDIR := .git
#
#   git_version
if [ -a $(GITDIR) ] ; \
then \
     GIT_VERSION := gitVersion.hpp ; \
fi;
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
#   Build scripts
SCRIPTS := scripts
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

.PHONY: all clean distclean install uninstall

all: $(HEADERS)/$(GIT_VERSION) $(BINDIR)/$(EXE) helpFile.pdf

$(HEADERS)/$(GIT_VERSION): .git
	bash $(SCRIPTS)/makeGitVersion.sh

$(BINDIR)/$(EXE): $(OBJS)
	mkdir -p $(BINDIR)
	$(CXX) $(LDFLAGS) $? -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

helpFile.pdf : db/helpFile.man
	bash $(SCRIPTS)/updateMan.sh

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(EXE)
	rm -f helpFile.pdf

install: $(BINDIR)/$(EXE)
	cp $(BINDIR)/$(EXE) $(INSTALL_DIR)/$(EXE)

uninstall:
	rm -fv $(INSTALL_DIR)/$(EXE)

distclean: clean
