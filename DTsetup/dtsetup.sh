#!/bin/bash

invalidOptionMessage="is an invalid option! Exiting..."
format="std"

#get arguments
while ! [[ -z "$1" ]] ; do
	case $1 in
		"-f" | "--format")
			shift
			format="$1"
		;;
		*)
			echo "$1" $invalidOptionMessage
			exit
		;;
	esac
	shift
done

#genearte a folder in parent dir with the name $newDirName and copy all DTA-filter
#files into $nweDirName using specified format
newDirName=$(date +"%y-%m-%d_%H%M%S")_dtarraySetup
echo "Creating dtarraySetup folder in parent dirrectory..."
mkdir $newDirName
echo "Copying DTA filter files to $newDirName..."
case $format in
	"std")
		#itterate through all dirs one level below parent dir looking for
		#DTA_filter files
		for D in ./* ; do
			#skip new dir just made
			if ! [ "${D/dtarraySetup}" = "$D" ] ; then
				continue;
			fi
			#copy each DTA filter file into dtarraySetup folder with the name of the dir
			#that it is currently stored in
			if [ -d "$D" ] ; then
				cd $D
				if [ -f DTASelect-filter.txt ] ; then
					echo -e "\tAdding $D..."
					cp DTASelect-filter.txt ../$newDirName/"$D".dtafilter
				fi
				cd ..
			fi
		done
	;;
	"subdir")
		#itterate through all dirs one level beliw parent dir looking for
		#DTA filter files
		for D in * ; do
			#skip new dir just made
			if ! [ "${D/dtarraySetup}" = "$D" ] ; then
				continue;
			fi
			if [ -d "$D" ] ; then
				cd $D
				#if DTA filter files exists in $D, copy file into a dir in $newDirName
				#with the same name as $D
				if [ -f DTASelect-filter.txt ] ; then
					echo -e "\tAdding $D..."
					mkdir ../$newDirName/"$D"/
					cp DTASelect-filter.txt ../$newDirName/"$D"/
				fi
				cd ..
			fi
		done
	;;
	*)
		echo $format" is not a valid output format! Exiting..."
	;;
esac

echo Done
