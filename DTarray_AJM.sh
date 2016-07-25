
#editable paramaters
scriptWD="$HOME/scripts/DTarray_AJM/"
recompileMessage='DTarray_AJM source code recompiled.'
invalidOptionMessage="is an invalid option! Exiting..."
numParamsInParamsFile=3
paramsFileName="dtarray_ajm.params"

#default paramaters
input="standard"
output="standard"
includeUnique="0"
wd=$(pwd)
recompile=false
sampleNamePrefx=""
rewriteParams=false
filesFound=false
keepParams=false

function writeParams {
echo 'outputFormat='$1 >> ./$paramsFileName
echo 'includeUnique='$2 >> ./$paramsFileName
echo 'sampleNamePrefix='$3 >> ./$paramsFileName
}

#get arguements
while ! [[ -z "$1" ]] ; do
    case $1 in
        "-i" | "--in")
            shift
            input="$1"
            ;;
        "-o" | "--out")
            shift
            output="$1"
            ;;
        "-u" | "--Unique")
            includeUnique="1"
            ;;
        "-d" | "--directory" )
			shift
            wd="$1"
            ;;
		"-r" | "--recompile")
			recompile=true
			;;
		"-p" | "--prefix")
			shift
			sampleNamePrefix="$1"
			;;
		"-rw")
			rewriteParams=true
			;;
		"-k" | "--keep")
			keepParams=true
			;;
        *)
            echo "$1" $invalidOptionMessage
            exit
        ;;
    esac
    shift
done

#summarize params for user
echo "input = "$input
echo "output = "$output
echo "includeUnique = "$includeUnique
echo "wd = "$wd
echo "recompile = "$recompile
echo "sampleNamePrefix = "$sampleNamePrefix
echo "rewriteParams = "$rewriteParams
echo "keepParams = "$keepParams

#compile source code if necissary
cd $scriptWD/
if ! [[ -a a.out ]] ; then
    recompileMessage='DTarray_AJM object file not found. Recompiling from source code...'
    recompile=true
fi
if $recompile ; then
    g++ utils.h dtafilter.h DTarray_AJM.cpp
    echo $recompileMessage
fi

#create params file
cd $wd
if $rewriteParams ; then
	mv $paramsFileName .dtarray_ajm.temp
fi
if ! [[ -a $paramsFileName ]] ; then
	echo -e "#Params for DTarray_AJM\n#File list generated on: "$(date +"%y-%m-%d_%H:%M:%S") > ./$paramsFileName
	echo -e "\n#Files to add" >> ./$paramsFileName
	case $input in
		"standard")
			#check if wd contains .dtafilter files
			if [ $(ls -l $wd/*.dtafilter | wc -l) -lt 1 ] ; then
				echo "DTA-filter files could not be found in the specified directory! Exiting..."
				exit
			fi
			#add all .dtafilter files to params file
			for f in *.dtafilter ; do
				colName=${f::${#f}-10}
				echo -e $colName'\t'$f >> $paramsFileName
			done
		;;
		"subdir")
			#check all directories on level below parent for DTASelect-filter files and
			#add all found to params file
			for D in * ; do
				if [ -d "$D" ] ; then
					cd "$D"
					if [ -a DTASelect-filter.txt ] ; then
						echo -e $D'\t'$D'/'"DTASelect-filter.txt" >> ../$paramsFileName
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
		;;
		*)
			echo "$input is not a valid input format! Exiting..."
			exit
		;;
	esac
	#add sampleNamePrefix to params
	echo -e "\n#Params" >> ./$paramsFileName
	writeParams $output $includeUnique $sampleNamePrefix
	fi
if [ -a $paramsFileName ] && ! $keepParams ; then
	echo "Params overwritten"
	numLines=$(($(echo $(wc -l $paramsFileName)|cut -d' ' -f 1)-$numParamsInParamsFile))
	head -n $numLines $paramsFileName > .dtarray_ajm.temp
	cat .dtarray_ajm.temp > $paramsFileName
	writeParams $output $includeUnique $sampleNamePrefix
fi

#run DTarray_AJM
cd $scriptWD
./a.out $wd/

echo "Done"
