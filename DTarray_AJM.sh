
#editable paramaters
recompile=false
scriptWD="$HOME/scripts/DTarray_AJM/"
recompileMessage='DTarray_AJM source code recompiled.'

#check if output format was specified. If not, use standard.
if [ -z "$1" ] ; then
    output="standard"
else output="$1"
fi

#check if dirrectory was specified. If not, use working dirrectory
if [ -z "$2" ] ; then
    wd=$(pwd)
else wd="$2"
fi

#compile source code if necissary
echo 'recompile = '$recompile
cd $scriptWD/
if ! [[ -a a.out ]] ; then
    recompileMessage='DTarray_AJM object file not found. Recompiling from source code...'
    recompile=true
fi
if $recompile ; then
    g++ DTarray_AJM.cpp
    echo $recompileMessage
fi

# If no params file is found in wd then create one by searching
# all dirrectories one level below wd for a DTASelect-filter file
# and add all folders where one is found to params file
cd $wd
if ! [[ -a dtarray_ajm.params ]] ; then
    for D in * ; do
        if [ -d "$D" ] ; then
            cd "$D"
            if [ -a DTASelect-filter.txt ] ; then
                echo -e $D'\t'$D'/'>> ../dtarray_ajm.params
            fi
        fi
        cd ..
    done
fi

#run DTarray_AJM
cd $scriptWD
./a.out $wd/ $output

echo "Done"