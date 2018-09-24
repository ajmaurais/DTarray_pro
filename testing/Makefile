
#commands
PYTHON := python

#dirs
TEMP_DIR := temp
DATA_DIR := data
SRC_DIR := src

#file names
SINGLE_ARG_LIST := test_arg_list.tsv
SINGLE_TEST_ARGS := test_args_single.txt
ALL_ARG_LIST := test_arg_list_short.tsv
ALL_TEST_ARGS := test_args_all.txt

.PHONY: argLists clean

argLists: $(TEMP_DIR)/$(SINGLE_TEST_ARGS) $(TEMP_DIR)/$(ALL_TEST_ARGS)

$(TEMP_DIR)/$(SINGLE_TEST_ARGS) : $(DATA_DIR)/$(SINGLE_ARG_LIST)
	mkdir -p $(TEMP_DIR)
	$(PYTHON) $(SRC_DIR)/makeArgList.py -c s $(DATA_DIR)/$(SINGLE_ARG_LIST) $(TEMP_DIR)/$(SINGLE_TEST_ARGS)

$(TEMP_DIR)/$(ALL_TEST_ARGS) : $(DATA_DIR)/$(ALL_ARG_LIST)
	mkdir -p $(TEMP_DIR)
	$(PYTHON) $(SRC_DIR)/makeArgList.py -c p $(DATA_DIR)/$(ALL_ARG_LIST) $(TEMP_DIR)/$(ALL_TEST_ARGS)

clean:
	rm -fv temp/*

