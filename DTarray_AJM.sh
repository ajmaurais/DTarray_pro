
#editable paramaters
scriptWD="$HOME/scripts/DTarray_AJM/"
recompileMessage='DTarray_AJM source code recompiled.'
invalidOptionMessage="is an invalid option! Exiting..."
numParamsInParamsFile=3
defaultParamsName="dtarray_ajm.params"
paramsCommentSymbol="#" #if changed, value must also be changed in utils.h

#default paramaters
paramsName=$defaultParamsName
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
	echo 'outputFormat='$1 >> ./$paramsName
	echo 'includeUnique='$2 >> ./$paramsName
	echo 'sampleNamePrefix='$3 >> ./$paramsName
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
		"-par")
			shift
			paramsName="$1"
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
echo -e "\nThe folowing paramaters were used:"
echo "input = "$input
echo "output = "$output
echo "includeUnique = "$includeUnique
echo "wd = "$wd
echo "recompile = "$recompile
echo "sampleNamePrefix = "$sampleNamePrefix
echo "rewriteParams = "$rewriteParams
echo "keepParams = "$keepParams
echo "paramsName = "$paramsName
echo

#compile source code if necissary
cd $scriptWD/
if ! [[ -a a.out ]] ; then
    recompileMessage='DTarray_AJM object file not found. Recompiling from source code...'
    recompile=true
fi
if $recompile ; then
	g++ main.cpp #>& compiler_output.txt
    echo $recompileMessage
fi

#create params file
cd $wd
if $rewriteParams ; then
	mv $paramsName .dtarray_ajm.temp
fi
if ! [[ -a $paramsName ]] ; then
	echo "Generating $paramsName using $input input format."
	echo -e "$paramsCommentSymbol Params for DTarray_AJM\n$paramsCommentSymbol File list generated on: "$(date +"%y-%m-%d_%H:%M:%S") > ./$paramsName
	echo -e "\n$paramsCommentSymbol Files to add" >> ./$paramsName
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
				echo -e $colName'\t'$f >> $paramsName
			done
		;;
		"subdir")
			#check all directories on level below parent for DTASelect-filter files and
			#add all found to params file
			for D in * ; do
				if [ -d "$D" ] ; then
					cd "$D"
					if [ -a DTASelect-filter.txt ] ; then
						echo -e $D'\t'$D'/'"DTASelect-filter.txt" >> ../$paramsName
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
	echo -e "\n$paramsCommentSymbol Params" >> ./$paramsName
	writeParams $output $includeUnique $sampleNamePrefix
	keepParams=true
fi
if [ -a $paramsName ] && ! $keepParams ; then
	numLines=$(($(echo $(wc -l $paramsName)|cut -d' ' -f 1)-$numParamsInParamsFile))
	head -n $numLines $paramsName > .dtarray_ajm.temp
	cat .dtarray_ajm.temp > $paramsName
	writeParams $output $includeUnique $sampleNamePrefix
fi

#run DTarray_AJM
cd $scriptWD
./a.out $wd/ $paramsName

echo "Done"
