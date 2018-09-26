
#
# C++ Compiler
CXX := g++
#
# latexmk
TEX := /Library/TeX/texbin/latexmk
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
#ifneq ("$(wildcard $(GITDIR))","")
GIT_VERSION := gitVersion.hpp
#endif
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
#   Tex dir
TEX_DIR := doc/tex
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

.PHONY: all gitVersion clean distclean install uninstall

all: gitVersion $(BINDIR)/$(EXE) helpFile.pdf DTarray_pro-Userguide.pdf $(BINDIR)/DTsetup

gitVersion :
	bash $(SCRIPTS)/makeGitVersion.sh

DTarray_pro-Userguide.pdf : $(TEX_DIR)/DTarray_pro-Userguide.tex
	cd $(TEX_DIR) && $(TEX) -pdf DTarray_pro-Userguide.tex
	cp $(TEX_DIR)/DTarray_pro-Userguide.pdf .

$(BINDIR)/$(EXE): $(OBJS)
	mkdir -p $(BINDIR)
	$(CXX) $(LDFLAGS) $? -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

helpFile.pdf : db/helpFile.man
	bash $(SCRIPTS)/updateMan.sh

$(BINDIR)/DTsetup : DTsetup/dtsetup.sh
	cp DTsetup/dtsetup.sh $(BINDIR)/DTsetup
	chmod +x $(BINDIR)/DTsetup

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(EXE) $(BINDIR)/DTsetup
	rm -f helpFile.pdf
	cd $(TEX_DIR) && rm -f ./*.aux ./*.dvi ./*.fdb_latexmk ./*.fls ./*.log ./*.out ./*.pdf ./*.toc 

install: $(BINDIR)/$(EXE)
	cp $(BINDIR)/$(EXE) $(INSTALL_DIR)/$(EXE)

uninstall:
	rm -fv $(INSTALL_DIR)/$(EXE)

distclean: clean
