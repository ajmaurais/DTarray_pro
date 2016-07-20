
#editable paramaters
scriptWD="$HOME/scripts/DTarray_AJM/"
recompileMessage='DTarray_AJM source code recompiled.'
invalidOptionMessage="is an invalid option! Exiting..."
numParamsInParamsFile=1

#default paramaters
input="standard"
output="standard"
includeUnique="0"
wd=$(pwd)
recompile=false
parseSamplePrefix="0"
sampleNamePrefix=""
rewriteParams=false
filesFound=false

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
        "-u" | "--uniuque")
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
			parseSamplePrefix="1"
			shift
			sampleNamePrefix="$1"
			;;
		"-rw")
			rewriteParams=true
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
echo "includeUniuque = "$includeUnique
echo "wd = "$wd
echo "recompile = "$recompile
echo "parseSamplePrefix = "$parseSamplePrefix
echo "sampleNamePrefix = "$sampleNamePrefix
echo "rewriteParams = "$rewriteParams

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
	mv dtarray_ajm.params .dtarray_ajm.temp
fi
if ! [[ -a dtarray_ajm.params ]] ; then
	echo -e "#Params for DTarray_AJM\n#File list generated on: "$(date +"%y-%m-%d_%H:%M:%S") > ./dtarray_ajm.params
	echo -e "\n#Files to add" >> ./dtarray_ajm.params
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
				echo -e $colName'\t'$f >> dtarray_ajm.params
			done
		;;
		"subdir")
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
		;;
		*)
			echo "$input is not a valid input format! Exiting..."
			exit
		;;
	esac
	#add sampleNamePrefix to params
	echo -e "\n#Params" >> ./dtarray_ajm.params
	echo "sampleNamePrefix="$sampleNamePrefix >> ./dtarray_ajm.params
fi
if [ -a dtarray_ajm.params ] ; then
	numLines=$(($(echo $(wc -l dtarray_ajm.params)|cut -d' ' -f 1)-$numParamsInParamsFile))
	head -n $numLines dtarray_ajm.params > .dtarray_ajm.temp
	cat .dtarray_ajm.temp > dtarray_ajm.params
	echo "sampleNamePrefix="$sampleNamePrefix >> dtarray_ajm.params
fi

#run DTarray_AJM
cd $scriptWD
./a.out $wd/ $output $includeUnique $parseSamplePrefix

echo "Done"
