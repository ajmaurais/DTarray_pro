/*
 DTarray_AJM reads in a specified number of dtaselect-filter files and writes protein, their molecular
 weights and spectral count data to OF_NAME in the working dirrectory.
 
 Written by Aaron Maurais
 25 May 2016
 */

#include <iostream>
#include <vector>
#include "dtafilter.h"

using namespace std;

//begin main
int main (int argc, char *argv[])
{
	//check paramaters
	string wd = string(argv[1]);
	assert(dirExists(wd));
	string outputFormat = "standard";
	if (argc >= 3)
		{
		outputFormat = argv[2];
		if(outputFormat != "standard" && outputFormat != "DB")
			{
			cout << outputFormat << " is not a valid output format. Exiting..." << endl;
			return 0;
			}
		}
	int includeUnique = 0;
	if (argc >= 4)
		{
		string includeUniqueStr = argv[3];
		if (includeUniqueStr != "0" && includeUniqueStr != "1")
			{
			cout << includeUniqueStr << " is not a valid arguement for includeUnique. Exiting..." << endl;
			return 0;
			}
		includeUnique = toInt(includeUniqueStr);
		}
	int parseSampleName = 0;
	if (argc >= 5)
		{
		string parseSampleNameStr = argv[4];
		if (parseSampleNameStr != "0" && parseSampleNameStr != "1")
			{
			cout << parseSampleNameStr << " is not a valid arguement for parseSampleName. Exiting..." << endl;
			return 0;
			}
		parseSampleName = toInt(parseSampleNameStr);
		}
	
	//read in names of files to combine from params file
	FilterFileParams filterFileParams;
	if (!filterFileParams.readDTParams(DTA_PARAMS_NAME, wd))
		{
		cout <<"Failed to read params file! Exiting..." << endl;
		return 0;
		}
	if(parseSampleName)
		cout << "Parsing colnames by prefix: " << filterFileParams.sampleNamePrefix << endl;
	
	//combine files
	cout << endl;
	Proteins proteins;
	proteins.initialize(filterFileParams);
	for (int i = 0; i < filterFileParams.numFiles; i++)
		{
		if(!proteins.readIn(wd+filterFileParams.file[i].path, filterFileParams.file[i].colname, includeUnique))
			{
			cout <<"Failed to read in " << filterFileParams.file[i].path <<"!" << endl <<
			"Exiting..." << endl;
			return 0;
			}
		cout << "Adding " << filterFileParams.file[i].colname << "..." << endl;
		}
	
	//write out combined data to OF_NAME
	if (outputFormat == "standard")
		{
		if(!proteins.writeOut(wd + OF_NAME, includeUnique, parseSampleName, filterFileParams.sampleNamePrefix))
			{
			cout << "Could not write outFile! Exiting..." << endl;
			return 0;
			}
		}
	else if (outputFormat == "DB")
		{
		if(!proteins.writeOutDB(wd + OF_NAME, includeUnique, parseSampleName, filterFileParams.sampleNamePrefix))
			{
			cout << "Could not write outFile! Exiting..." << endl;
			return 0;
			}
		cout << endl << "Results written in database format." << endl;
		}
	else
		{
		cout << endl << outputFormat << " is not a valid output format! Exiting..." << endl;
		return 0;
		}
	
	//summarize results for user
	cout << proteins.colNames.size() << " files combined." << endl;
	cout << "Results written to: " << OF_NAME << endl;
	
	return 0;
}
//end main

