
#editable paramaters
recompile=false
scriptWD="$HOME/scripts/DTarray_AJM/"
recompileMessage='DTarray_AJM source code recompiled.'

#check if output format was specified. If not, use standard.
if [ -z "$1" ] ; then
    output="standard"
else output="$1"
fi

#check if user wants to include unique peptide spectral counts
#by default answer is no
if [ -z "$2" ] ; then
    includeUnique="0"
else includeUnique="$2"
fi

#check if dirrectory was specified. If not, use working dirrectory
if [ -z "$3" ] ; then
    wd=$(pwd)
else wd="$3"
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

#check if wd contains .ms2 files
if [ $(ls -lR $wd/*.dtafilter | wc -l) -lt 1 ] ; then
    echo "DTA-filter files could not be found in the specified directory! Exiting..."
    exit
fi

# If no params file is found in wd then create one by searching
# all dirrectories one level below wd for a DTASelect-filter file
# and add all folders where one is found to params file
cd $wd
if ! [[ -a dtarray_ajm.params ]] ; then
    for f in *.dtafilter ; do
        colName=${f::${#f}-10}
        echo -e $colName'\t'$f >> dtarray_ajm.params
    done
fi

#run DTarray_AJM
cd $scriptWD
./a.out $wd/ $output $includeUnique

echo "Done"