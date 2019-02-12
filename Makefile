
#
# C++ Compiler
CXX := g++
#
# latexmk
TEX := $(shell command -v latexmk 2> /dev/null)
#
# git version vars
GITVERSION:= $(shell git log -1 --pretty='%h')
GITDATE:= $(shell git log -1 --format=%cd --date=local)
#
# Flags
#
#   Compiler
CXXFLAGS += -c -g -Wall -std=c++11 -DGIT_VERSION="\"${GITVERSION}\"" -DGIT_DATE="\"${GITDATE}\""
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
ifneq ("$(wildcard $(GITDIR))","")
GIT_EXISTS = 1
else
GIT_EXISTS = 0
endif
#
#   git_version
GIT_VERSION := gitVersion.hpp
#
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

.PHONY: all clean distclean install uninstall

#TARGETS = $(HEADERDIR)/$(GIT_VERSION) $(BINDIR)/$(EXE) $(BINDIR)/DTsetup helpFile.pdf DTarray_pro-Userguide.pdf
TARGETS = $(BINDIR)/$(EXE) $(BINDIR)/DTsetup helpFile.pdf DTarray_pro-Userguide.pdf

all: $(TARGETS)

DTarray_pro-Userguide.pdf : $(TEX_DIR)/DTarray_pro-Userguide.tex
ifndef TEX
	$(warning "No latexmk in $(PATH), skipping build of DTarray_pro-Userguide.pdf")
else
	cd $(TEX_DIR) && latexmk -pdf DTarray_pro-Userguide.tex
	cp $(TEX_DIR)/DTarray_pro-Userguide.pdf .
endif

$(BINDIR)/$(EXE): $(OBJS)
	mkdir -p $(BINDIR)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERDIR)/%.hpp
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
