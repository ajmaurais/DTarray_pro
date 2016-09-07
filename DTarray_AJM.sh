
#editable paramaters
scriptWDHome="$HOME/scripts/DTarray_AJM"
scriptWDbin="$scriptWDHome/bin"
scriptWDsrc="$scriptWDHome/src"
scriptWDdb="$scriptWDHome/db"
locDBfname="$scriptWDdb/subCelluarLoc.txt"
aaDBfanme="$scriptWDdb/aaMasses.txt"
staticModificationsDB="$scriptWDdb/staticModifications.txt"
binName="DTarray_AJM"

recompileMessage='DTarray_AJM source code recompiled.'
invalidOptionMessage="is an invalid option! Exiting..."
defaultParamsName="dtarray_ajm.params"
defaultFlistName="dtarray_ajm_flist.txt"
defaultStaticModificationsName="staticModifications.txt"
paramsCommentSymbol="#" #if changed, COMMENT_SYMBOL must also be changed in src/utils.h

#default paramaters
paramsName=$defaultParamsName
flistName=$defaultFlistName
staticModificationsName=$defaultStaticModificationsName
input="standard"
output="standard"
includeUnique="0"
wd=$(pwd)'/'
recompile=false
sampleNamePrefx=""
rewriteFlist=false
filesFound=false
keepParams=false
getSubCelluarLoc="0"
calcMW=false
calcMWStr="1"
peptideDBfname=""
rewriteSmod=false

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
			shift
			arg="$1"
			case $arg in
				"flist")
					rewriteFlist=true
				;;
				"smod")
					rewriteSmod=true
				;;
			esac
			;;
		"-k" | "--keep")
			keepParams=true
			;;
		"-par")
			shift
			flistName="$1"
			keepParams=true
			;;
		"-loc")
			getSubCelluarLoc="1"
			;;
		"-mw")
			shift
			peptideDBfname="$1"
			calcMWStr="1"
			calcMW=true
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
echo "keepParams = "$keepParams
echo "flistName = "$flistName
echo "paramsName = "$paramsName
echo "getSubCelluarLoc = "$getSubCelluarLoc
echo "calcMW = "$calcMW
echo "peptideDBfname = "$peptideDBfname
echo "rewriteSmod = "$rewriteSmod
echo

#compile source code if necissary
cd $scriptWDbin
if ! [[ -a $binName ]] ; then
    recompileMessage='DTarray_AJM object file not found. Recompiling from source code...'
    recompile=true
fi
if $recompile ; then
	cd $scriptWDsrc
	g++ -o $scriptWDbin/$binName main.cpp
    echo $recompileMessage
fi

#create params file
cd $wd
if $rewriteFlist ; then
	mv $flistName .dtarray_ajm.temp
fi
if ! [[ -a $flistName ]] ; then
	echo "Generating $flistName using $input input format."
	echo -e "$paramsCommentSymbol File list for DTarray_AJM\n$paramsCommentSymbol File list generated on: "$(date +"%y-%m-%d_%H:%M:%S")'\n' > ./$flistName
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
				echo -e $colName'\t'$f >> $flistName
			done
		;;
		"subdir")
			#check all directories on level below parent for DTASelect-filter files and
			#add all found to params file
			for D in * ; do
				if [ -d "$D" ] ; then
					cd "$D"
					if [ -a DTASelect-filter.txt ] ; then
						echo -e $D'\t'$D'/'"DTASelect-filter.txt" >> ../$flistName
						filesFound=true
					fi
					cd ..
				fi
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
fi

#write output paramaters to params file
cd $wd
if ! $keepParams ; then
	echo -e "$paramsCommentSymbol Params for DTarray_AJM\n$paramsCommentSymbol Params generated on: "$(date +"%y-%m-%d_%H:%M:%S")'\n' > ./$paramsName
	echo -e '<paramsFile>\n' >> ./$paramsName
	echo -e '<params>\n' >> ./$paramsName
	echo 'outputFormat='$output >> ./$paramsName
	echo 'includeUnique='$includeUnique >> ./$paramsName
	echo 'sampleNamePrefix='$sampleNamePrefix >> ./$paramsName
	echo 'getSubCelluarLoc='$getSubCelluarLoc >> ./$paramsName
	echo 'locDBfname='$locDBfname >> ./$paramsName
	echo 'calcMW='$calcMWStr >> ./$paramsName
	echo 'peptideDBfname='$peptideDBfname >> ./$paramsName
	echo 'aaDBfanme='$aaDBfanme >> ./$paramsName
	echo 'staticModsFname='$staticModificationsName >> ./$paramsName
	echo -e '\n</params>\n' >> ./$paramsName
	echo -e '</paramsFile>' >> ./$paramsName
fi

if $calcMW == "1" ; then
	if ! [ -a $staticModificationsName ] || $rewriteSmod ; then
		echo -e "$paramsCommentSymbol Static modifications for DTarray_AJM\n$paramsCommentSymbol file generated on: "$(date +"%y-%m-%d_%H:%M:%S")'\n' > $staticModificationsName
		echo -e '<staticModifications>\n' >> ./$staticModificationsName
		cat $staticModificationsDB >> ./$staticModificationsName
		echo -e '\n\n</staticModifications>\n' >> ./$staticModificationsName
	fi
fi


#run DTarray_AJM
cd $scriptWDbin
#./$binName $wd $flistName $paramsName

echo "Done"
