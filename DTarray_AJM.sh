
#editable paramaters
recompile=false
scriptWD="$HOME/scripts/DTarray_AJM/"
recompileMessage='DTarray_AJM source code recompiled.'
invalidOptionMessage="is an invalid option! Exiting..."

#default paramaters
input="standard"
output="standard"
includeUnique="0"
wd=$(pwd)
filesFound=false

#get arguements
while ! [[ -z "$1" ]] ; do
    case $1 in
        "-in")
            shift
            input="$1"
            ;;
        "-out")
            shift
            output="$1"
            ;;
        "-uniuque")
            includeUnique="1"
            ;;
        "-dir")
            shift
            wd="$1"
            ;;
        *)
            echo "$1" $invalidOptionMessage
            exit
        ;;
    esac
    shift
done

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

#create params file
cd $wd
case $input in
    "standard")
        if ! [[ -a dtarray_ajm.params ]] ; then
            #check if wd contains .dtafilter files
            if [ $(ls -l $wd/*.dtafilter | wc -l) -lt 1 ] ; then
                echo "DTA-filter files could not be found in the specified directory! Exiting..."
                exit
            fi
			#add all .dtafilter files to params file
            for f in *.dtafilter ; do
                colName=${f::${#f}-10}
                echo -e $colName'\t'$f >> dtarray_ajm.params
            done
        fi
	;;
    "subdir")
        if ! [[ -a dtarray_ajm.params ]] ; then
			#check all directories on level below parent for DTASelect-filter files and
			#add all found to params file
			for D in * ; do
                if [ -d "$D" ] ; then
                    cd "$D"
                    if [ -a DTASelect-filter.txt ] ; then
                        echo -e $D'\t'$D'/'"DTASelect-filter.txt" >> ../dtarray_ajm.params
                        filesFound=true
                    fi
                fi
                cd ..
            done
			#if no DTA-filter files were found, exit program.
			if ! $filesFound ; then
				echo "No DTA-filter files were found! Exiting..."
				exit
			fi
		fi
	;;
	*)
		echo "$input is not a valid input format! Exiting..."
		exit
	;;
esac


#run DTarray_AJM
cd $scriptWD
./a.out $wd/ $output $includeUnique

echo "Done"
