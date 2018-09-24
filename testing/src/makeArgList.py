# -*- coding: utf-8 -*-

import itertools
import csv
import argparse
import sys

COLS_OPTION = 'option'
COLS_ARG = 'arg'
REQUIRED_COLS = [COLS_OPTION, COLS_ARG]

def checkCols(_colNames):
    for col in REQUIRED_COLS:
        if col not in _colNames.keys():
            raise RuntimeError('{} is a required column!'.format(col))

def main(argv):

    parser = argparse.ArgumentParser(prog = 'makeArgList')
    
    parser.add_argument('arg_list',
                        help = '.tsv file with all arguments to consider.')
    
    parser.add_argument('output', help = 'Output file name')
    
    parser.add_argument('-a', '--add_dash',
                        type = int,
                        choices = [0, 1],
                        default = 1)
    
    args = parser.parse_args()
    
    #read arg_list
    with open(args.arg_list, 'r') as f:
        reader = csv.reader(f)
        lines = list(reader)
    
    #parse arg_list
    argLists = dict()
    colNames = dict()
    for i, line in enumerate(lines):
        if i == 0:
            colNames = {col:i for i, col in enumerate(line[0].split('\t'))}
            checkCols(colNames)
        else:
            cols = line[0].split('\t')
            
            arg = cols[colNames[COLS_ARG]].strip()
            if args.add_dash and arg != '':
                arg = '-' + arg
            
            if cols[colNames[COLS_OPTION]] in argLists.keys():
                argLists[cols[colNames[COLS_OPTION]]].append(arg)
            else:
                argLists[cols[colNames[COLS_OPTION]]] = [arg]
    
    prod = 0
    for key in argLists.keys():
        length = len(argLists[key])
        print('{} : {}'.format(key, length))
        if prod == 0:
            prod = length
        else: prod *= length
    
    print('product is {}'.format(prod))
    
    #create list of possible combinations of args
    tests = list()
    for comb in itertools.product(*argLists.values()):        
        #conc args
        arg = ''
        for x in comb:
            if x == '':
                continue
            arg = arg + ' ' + x

        tests.append(arg)
        
    outF = open(args.output, 'w')
    
    for arg in tests:
        outF.write('{}\n'.format(arg))
    
if __name__ == '__main__':
    main(sys.argv)
