
#commands
PYTHON := python3

#dirs
SRC_DIR := src
DTARRAY_BIN := ../bin/DTarray

#file names
SINGLE_ARG_LIST := test_arg_list_single.tsv
SINGLE_TEST_ARGS := test_args_single.txt
ALL_ARG_LIST := test_arg_list_short.tsv
ALL_TEST_ARGS := test_args_all.txt
TEST_SCRIPT := $(SRC_DIR)/runTests.py

.PHONY: argLists errorTest clean

argLists: temp/$(SINGLE_TEST_ARGS) temp/$(ALL_TEST_ARGS)

errorTest : temp/$(SINGLE_TEST_ARGS)
	$(PYTHON) $(TEST_SCRIPT) -v temp/test_args_single.txt $(DTARRAY_BIN)

temp/$(SINGLE_TEST_ARGS) : testLists/$(SINGLE_ARG_LIST)
	mkdir -p temp
	$(PYTHON) $(SRC_DIR)/makeArgList.py -c s testLists/$(SINGLE_ARG_LIST) temp/$(SINGLE_TEST_ARGS)

temp/$(ALL_TEST_ARGS) : testLists/$(ALL_ARG_LIST)
	mkdir -p temp
	$(PYTHON) $(SRC_DIR)/makeArgList.py -c p testLists/$(ALL_ARG_LIST) temp/$(ALL_TEST_ARGS)



clean:
	rm -fv temp/*
	rm -rfv log/*
