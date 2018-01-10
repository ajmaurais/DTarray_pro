#/bin/bash

rsync -vr -e 'ssh -p 22022' --exclude='.*' build.sh db include src mauraisa@pleiades.bc.edu:~/scripts/DTarray_pro/
