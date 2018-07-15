#/bin/bash

rsync -vr -e 'ssh -p 22022' --exclude='.*' Makefile helpFile.pdf db include src mauraisa@pleiades.bc.edu:~/local/DTarray_pro/
