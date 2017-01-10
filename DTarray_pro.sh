
#editable paramaters
scriptWDHome="$HOME/scripts/DTarray_pro"
scriptWDbin="$scriptWDHome/bin"
scriptWDsrc="$scriptWDHome/src"
scriptWDdb="$scriptWDHome/db"
locDBfname="$scriptWDdb/subCelluarLoc.txt"
seqDBfname="$scriptWDdb/humanProteome.fasta"
aaDBfanme="$scriptWDdb/aaMasses.txt"
fxnDBfname="$scriptWDdb/humanFxn.txt"
staticModificationsDB="$scriptWDdb/staticModifications.txt"
binName="DTarray_pro"
helpFileFname="$scriptWDdb/helpFile.man"
versionNum='1.4'

recompileMessage='DTarray_pro source code recompiled.'
invalidOptionMessage="is an invalid option! Exiting...\nUse DTarray -h for help."
defaultParamsName="dtarray_pro.params"
defaultFlistName="dtarray_pro_flist.txt"
defaultStaticModificationsName="staticModifications.txt"
paramsCommentSymbol="#" #if changed, COMMENT_SYMBOL must also be changed in src/DTarray_AJM.hpp

#default paramaters
paramsName=$defaultParamsName
flistName=$defaultFlistName
staticModificationsName=$defaultStaticModificationsName
ofname="DTarray_pro.tsv"
dbOfname="DTarray_long.tsv"
peptideOfFname="peptideList.tsv"
dbPeptideOfFname="peptideList_long.tsv"
input="standard"
output="1"
includeUnique="0"
wd=$(pwd)'/'
recompile=false
sampleNamePrefx=""
rewriteFlist=false
filesFound=false
keepParams=false
getSubCelluarLoc="0"
seqDBFnameTr=""
calcMW=false
calcMWStr="0"
mwDBFname=""
rewriteSmod=false
getSeq="0"
continueExc=true
includePeptides="0"
peptideOutput="0"
includeCoverage="0"
getFxn="0"
useDefaultSeqDB="1"
includeNullPeptides="0"
supInfoOutput="0"
groupPeptides="1"

function usage {
	cat $scriptWDdb/usage.txt
	echo
	exit
}

function isArg {
	if [[ $1 == -* ]] || [[ -z "$1" ]] ; then
		usage
	fi
}

function absPath {
	if [ -d "$1" ]; then
		( cd "$1"; dirs -l +0 )
	else
		( cd "$(dirname "$1")"; d=$(dirs -l +0); echo "${d%/}/${1##*/}" )
	fi
}

function purgeDir {
	rm -v ./$defaultParamsName
	rm -v ./$defaultFlistName
	rm -v ./$staticModificationsName
	rm -v ./$ofname
	rm -v ./$dbOfname
	rm -v ./$peptideOfFname
	rm -v ./$dbPeptideOfFname
}

#get arguements
while ! [[ -z "$1" ]] ; do
    case $1 in
		"-h" | "--help")
			man $helpFileFname
			exit;;
		"-of" | "--ofname")
			shift
			isArg "$1"
			ofname="$1";;
        "-i" | "--in")
            shift
			isArg "$1"
            input="$1";;
        "-o" | "--out")
            shift
			isArg "$1"
            output="$1";;
        "-u" | "--unique")
            includeUnique="1";;
        "-d" | "--dir" )
			shift
			isArg "$1"
            wd=$(absPath "$1")
			echo $wd;;
		"-r" | "--recompile")
			recompile=true
			continueExc=false;;
		"-f" | "--prefix")
			shift
			isArg "$1"
			sampleNamePrefix="$1";;
		"-rw")
			shift
			isArg "$1"
			arg="$1"
			case $arg in
				"flist")
					rewriteFlist=true;;
				"smod")
					rewriteSmod=true;;
				*)
					echo -e "$1" $invalidOptionMessage
					exit;;
			esac
			;;
		"-k" | "--keepParams")
			keepParams=true;;
		"-par")
			shift
			isArg "$1"
			flistName="$1"
			keepParams=true;;
		"-loc")
			getSubCelluarLoc="1";;
		"-mw")
			if [[ $2 == -* ]] || [[ -z "$2" ]] ; then
				mwDBFname=$seqDBfname
				calcMWStr="1"
				calcMW=true
				useDefaultSeqDB="0"
			else
				shift
				isArg '$1'
				mwDBFname=$(absPath "$1")
				calcMWStr="1"
				calcMW=true
				useDefaultSeqDB="1"
			fi;;
		"-seq")
				if [[ $2 == -* ]] || [[ -z "$2" ]] ; then
					seqDBFnameTr=$seqDBfname
					getSeq="1"
				else
					shift
					isArg "$1"
					seqDBFnameTr=$(absPath "$1")
					getSeq="1"
				fi;;
		"-p" | "--peptides")
			shift
			isArg "$1"
			peptideOutput="$1"
			includePeptides=true
			;;
		"-c" | "--coverage")
			includeCoverage="1" ;;
		"-s")
			shift
			isArg "$1"
			supInfoOutput="$1" ;;
		"--purge")
			purgeDir
			exit ;;
		"-pswd")
			echo $scriptWDHome
			exit ;;
		"-cdswd")
			cd $scriptWDHome
			exit ;;
		"-oswd")
			open $scriptWDHome
			exit ;;
		"-fxn")
			getFxn="1" ;;
		"-g" | "--group")
			shift
			isArg "$1"
			groupPeptides="$1" ;;
        *)
            echo -e "$1" $invalidOptionMessage
			usage;;
    esac
    shift
done

#compile source code if necissary
cd $scriptWDbin
if ! [[ -a $binName ]] ; then
    recompileMessage='DTarray_AJM object file not found. Recompiling from source code...'
    recompile=true
	continueExc=true
fi
if $recompile ; then
	cd $scriptWDsrc
	g++ -o $scriptWDbin/$binName main.cpp
    echo $recompileMessage
	if ! $continueExc ; then
		exit
	fi
fi

echo -e '\nDTarray_pro v'$versionNum
echo

#create params file
cd $wd
if $rewriteFlist ; then
	mv $flistName .dtarray_ajm.temp
fi
if ! [[ -a $flistName ]] ; then
	echo "Generating $flistName using $input input format."
	echo -e "$paramsCommentSymbol File list for DTarray_pro\n$paramsCommentSymbol File list generated on: "$(date +"%y-%m-%d_%H:%M:%S")'\n' > ./$flistName
	echo -e "\n<versionNum>"$versionNum"</versionNum>\n<flist>\n" >> ./$flistName
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
	echo -e "\n</flist>\n" >> ./$flistName
fi

#write output paramaters to params file
cd $wd
if ! $keepParams ; then
	echo -e "$paramsCommentSymbol Params for DTarray_pro\n$paramsCommentSymbol Params generated on: "$(date +"%y-%m-%d_%H:%M:%S")'\n' > ./$paramsName
	echo -e '<paramsFile>\n' >> ./$paramsName
	echo -e '<params>\n' >> ./$paramsName
	echo 'versionNum'=$versionNum >> ./$paramsName
	echo 'ofname='$ofname >> ./$paramsName
	echo 'dbOfname='$dbOfname >> ./$paramsName
	echo 'peptideOfFname='$peptideOfFname >> ./$paramsName
	echo 'dbPeptideOfFname='$dbPeptideOfFname >> ./$paramsName
	echo 'outputFormat='$output >> ./$paramsName
	echo 'includeUnique='$includeUnique >> ./$paramsName
	echo 'sampleNamePrefix='$sampleNamePrefix >> ./$paramsName
	echo 'getSubCelluarLoc='$getSubCelluarLoc >> ./$paramsName
	echo 'locDBfname='$locDBfname >> ./$paramsName
	echo 'calcMW='$calcMWStr >> ./$paramsName
	echo 'mwDBFname='$mwDBFname >> ./$paramsName
	echo 'aaDBfanme='$aaDBfanme >> ./$paramsName
	echo 'getSeq='$getSeq >> ./$paramsName
	echo 'seqDBfname='$seqDBFnameTr >> ./$paramsName
	echo 'getFxn='$getFxn >> ./$paramsName
	echo 'fxnDBfname='$fxnDBfname >> ./$paramsName
	echo 'staticModsFname='$staticModificationsName >> ./$paramsName
	echo 'includePeptides='$peptideOutput >> ./$paramsName
	echo 'includeCoverage='$includeCoverage >> ./$paramsName
	echo 'useDefaultSeqDB='$useDefaultSeqDB >> ./$paramsName
	echo 'includeNullPeptides='$includeNullPeptides >> ./$paramsName
	echo 'supInfoOutput='$supInfoOutput >> ./$paramsName
	echo 'groupPeptides='$groupPeptides >> ./$paramsName
	echo -e '\n</params>\n' >> ./$paramsName
	echo -e '</paramsFile>' >> ./$paramsName
fi

if $calcMW ; then
	if ! [ -a $staticModificationsName ] || $rewriteSmod ; then
		echo -e "$paramsCommentSymbol Static modifications for DTarray_pro\n$paramsCommentSymbol file generated on: "$(date +"%y-%m-%d_%H:%M:%S")'\n' > $staticModificationsName
		echo -e '<staticModifications>\n' >> ./$staticModificationsName
		cat $staticModificationsDB >> ./$staticModificationsName
		echo -e '\n\n</staticModifications>\n' >> ./$staticModificationsName
	fi
fi

#run DTarray_pro
cd $scriptWDbin
./$binName $wd'/' $flistName $paramsName

echo "Done"
