# -*- coding: utf-8 -*-

import argparse
import sys
import subprocess as sp
import os
import datetime

def readArgList(fname, binary):
    inF = open(fname, 'r')
    lines = inF.readlines()
    
    for i in range(len(lines)):
        lines[i] = '{} {}'.format(binary, lines[i].strip())
    
    return lines

def setupLogDir(logPath):
    timeStamp = datetime.datetime.now().strftime('%y-%m-%d_%H%M%S')
    dirName = '{}/{}_log'.format(logPath, timeStamp)
    os.makedirs(dirName, exist_ok = True)
    return timeStamp, dirName

def writeSummary(out, binPath, testType, timeStamp, nTests, nError, printHeader = True):
    
    #get absolute path of binTested
    #first check if binary file exists in current folder
    if os.path.exists(os.path.abspath(binPath)):
        fullBinPath = os.path.abspath(binPath)
    else: #otherwise try to find binPath on PATH env variable
       p = sp.Popen(['which', binPath], stdout=sp.PIPE)
       stdout, stderr = p.communicate()
       assert(p.wait() == 0)
       fullBinPath = stdout.decode('utf-8').strip()
    
    #write file
    if printHeader:
        out.write('paramater\tvalue\n')
    
    out.write('paramater\tvalue\n'
               'binPath\t{}\n'
               'testType\t{}\n'
               'timeStamp\t{}\n'
               'nTests\t{}\n'
               'nError\t{}\n'.format(fullBinPath, testType, timeStamp, nTests, nError))
    

def main(argv):
    
    parser = argparse.ArgumentParser(prog = 'runTests',
                                     description = 'Run tests for DTaray')
    
    parser.add_argument('argList',
                        help = 'arg list generated by makeArgList to test.')
    
    parser.add_argument('binaryPath',
                        help = 'Path to binary to test')
    
    parser.add_argument('-l', '--logPath',
                        help = 'Path to print log files. Default is ./log',
                        default = 'log/')
    
    parser.add_argument('-d', '--dataPath',
                        help = 'Path to data folder. Default is ./data',
                        default = 'data/')
    
    parser.add_argument('-t', '--type',
                        choices = ['error', 'output'],
                        default = 'error',
                        help = 'Choose test type. error will check for runtime'
                        'errors. output will check outputs against validated'
                        'output files in testing/standards/. error is the default.')
    
    parser.add_argument('-p', '--progressStop',
                        help = 'Choose how often progress is updated.',
                        type = int,
                        default = 100)
    
    parser.add_argument('-v', '--verbose',
                        action = 'store_true',
                        default = False)
    
    args = parser.parse_args()
    
    #setup log files
    timeStamp, dirName = setupLogDir(args.logPath)
    all_outF = open('{}/all.txt'.format(dirName), 'w')
    error_outF = open('{}/error.txt'.format(dirName), 'w')
    
    #run tests
    tests = readArgList(args.argList, args.binaryPath)
    nTests = len(tests)
    nError = 0
    for i, test in enumerate(tests):
        
        if i % args.progressStop == 0:
            sys.stdout.write('Working on test {} of {}\n'.format(i, nTests))
        
        testCommand = '\nTesting: {}'.format(test)
        all_outF.write(testCommand)
        
        #call subprocess and get data
        p = sp.Popen(test, stdout=sp.PIPE, stderr = sp.PIPE, 
                     cwd = args.dataPath, shell = True)
        stdout, stderr = p.communicate()
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        pStatus = p.wait()
    
        if pStatus != 0:
            nError += 1
            error_outF.write(testCommand)
            lines = stderr.split('\n')
            for line in lines:
                if line:
                    l = '\n\t{}'.format(line)
                    all_outF.write(l)
                    error_outF.write(l)
                    
        else:
            all_outF.write(' -> Done')
    
    #write summary
    outF = open('{}/summary.txt'.format(dirName), 'w')
    writeSummary(outF, args.binaryPath, args.type, timeStamp, nTests, nError)
    
    if args.verbose:
        writeSummary(sys.stdout, args.binaryPath, args.type, timeStamp, nTests, nError)
    
    sys.stdout.write('Done\n')
    
if __name__ == '__main__':
    main(sys.argv)
