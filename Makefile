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
INCLUDEFLAGS := -I./utils/include
#
#
# Program name
EXE := DTarray
#
#
# Directories
#
# utils
UTILS_DIR := utils
UTILS_LIB := lib/utils.a
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
#	Libraries
LIBDIR := lib
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
#
################################################################################

HEADERS := $(wildcard $(HEADERDIR)/*.hpp)
SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(subst $(SRCDIR)/,$(OBJDIR)/,$(SRCS:.cpp=.o))

CXXFLAGS += $(INCLUDEFLAGS) -I$(HEADERDIR)
LDFLAGS += $(LIBFLAGS) -L$(LIBDIR)

.PHONY: all clean distclean doc

TARGETS = $(BINDIR)/$(EXE) $(BINDIR)/DTsetup helpFile.pdf DTarray_pro-Userguide.pdf

all: $(TARGETS)

DTarray_pro-Userguide.pdf : $(TEX_DIR)/DTarray_pro-Userguide.tex
ifndef TEX
	$(warning "No latexmk in $(PATH), skipping build of DTarray_pro-Userguide.pdf")
else
	cd $(TEX_DIR) && latexmk -pdf DTarray_pro-Userguide.tex
	cp $(TEX_DIR)/DTarray_pro-Userguide.pdf .
endif

$(BINDIR)/$(EXE): $(UTILS_LIB) $(OBJS)
	mkdir -p $(BINDIR)
	$(CXX) $(LDFLAGS) $(OBJS) $(UTILS_LIB) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(UTILS_LIB):
	cd $(UTILS_DIR); $(MAKE)
	mkdir -p $(LIBDIR)
	cd $(LIBDIR); ln -s ../$(UTILS_DIR)/$(UTILS_LIB)

helpFile.pdf : db/helpFile.man
	bash $(SCRIPTS)/updateMan.sh

$(BINDIR)/DTsetup : DTsetup/dtsetup.sh
	cp DTsetup/dtsetup.sh $(BINDIR)/DTsetup
	chmod +x $(BINDIR)/DTsetup

doc:
	doxygen doc/doxygen/Doxyfile

clean:
	rm -f $(OBJDIR)/*.o $(BINDIR)/$(EXE) $(BINDIR)/DTsetup $(LIBDIR)/*.a
	rm -f helpFile.pdf
	cd $(TEX_DIR) && rm -f ./*.aux ./*.dvi ./*.fdb_latexmk ./*.fls ./*.log ./*.out ./*.pdf ./*.toc
	cd $(UTILS_DIR); $(MAKE) clean

distclean: clean
